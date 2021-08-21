.. _cartopy: https://scitools.org.uk/cartopy/docs/latest/

.. _basemap: https://matplotlib.org/basemap/index.html

.. _seaborn: https://seaborn.pydata.org

.. _pandas: https://pandas.pydata.org

.. _xarray: http://xarray.pydata.org/en/stable/

.. _rainbow: https://doi.org/10.1175/BAMS-D-13-00155.1

.. _xkcd: https://blog.xkcd.com/2010/05/03/color-survey-results/

.. _opencolor: https://yeun.github.io/open-color/

.. _cmocean: https://matplotlib.org/cmocean/

.. _fabio: http://www.fabiocrameri.ch/colourmaps.php

.. _brewer: http://colorbrewer2.org/

.. _sciviscolor: https://sciviscolor.org/home/colormoves/

.. _matplotlib: https://matplotlib.org/stable/tutorials/colors/colormaps.html

.. _seacolor: https://seaborn.pydata.org/tutorial/color_palettes.html

.. _texgyre: https://frommindtotype.wordpress.com/2018/04/23/the-tex-gyre-font-family/

.. _why:

============
Why ProPlot?
============

Matplotlib is an extremely powerful plotting package used by
scientists and engineers far and wide. However,
matplotlib can be cumbersome or repetitive for users who...

* Make highly complex figures with many subplots.
* Want to finely tune their annotations and aesthetics.
* Need to make new figures nearly every day.

ProPlot's core mission is to provide a smoother plotting experience for
matplotlib's most demanding users. We accomplish this by *expanding upon*
matplotlib's :ref:`object-oriented interface <usage_background>`. ProPlot
makes changes that would be hard to justify or difficult to incorporate
into matplotlib itself, owing to differing design choices and backwards
compatibility considerations.

This page enumerates these changes and explains how they address the
limitations of matplotlib's default interface. To start using these
features, see the :ref:`usage introduction <usage>`
and the :ref:`user guide <ug_basics>`.

.. _why_less_typing:

Less typing, more plotting
==========================

.. rubric:: Limitation

Matplotlib users often need to change lots of plot settings all at once. With
the default interface, this requires calling a series of one-liner setter methods.

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
   interface to the object-oriented interface (see :ref:`Using ProPlot`).

.. rubric:: Solution

ProPlot introduces the `proplot.axes.Axes.format` command to resolve this.
Think of this as an expanded and thoroughly documented version of the
`matplotlib.artist.Artist.update` command. `~proplot.axes.Axes.format` can modify things
like axis labels and titles and apply new :ref:`"rc" settings <why_rc>` to existing
axes. It also integrates with various :ref:`constructor functions <why_constructor>`
to help keep things succinct. Further, the `proplot.figure.Figure.format`
and `proplot.figure.SubplotGrid.format` commands can be used to
`~proplot.axes.Axes.format` several subplots at once.

Together, these features significantly reduce the amount of code needed to create
highly customized figures. As an example, it is trivial to see that...

.. code-block:: python

   import proplot as pplt
   fig, axs = pplt.subplots(ncols=2)
   axs.format(linewidth=1, color='gray')
   axs.format(xlim=(0, 100), xticks=10, xtickminor=True, xlabel='foo', ylabel='bar')

is much more succinct than...

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import matplotlib as mpl
   with mpl.rc_context(rc={'axes.linewidth': 1, 'axes.color': 'gray'}):
       fig, axs = plt.subplots(ncols=2, sharey=True)
       axs[0].set_ylabel('bar', color='gray')
       for ax in axs:
           ax.set_xlim(0, 100)
           ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
           ax.tick_params(width=1, color='gray', labelcolor='gray')
           ax.tick_params(axis='x', which='minor', bottom=True)
           ax.set_xlabel('foo', color='gray')

.. _why_constructor:

Class constructor functions
===========================

.. rubric:: Limitation

Matplotlib and `cartopy`_ define several classes with verbose names like
`~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`,
and `~cartopy.crs.LambertAzimuthalEqualArea`. They also keep them out of the top-level
package namespace. Since plotting code has a half life of about 30 seconds, typing out
these extra class names and import statements can be a *major* drag.

Parts of matplotlib's interface were actually designed with this in mind.
`Backend classes <https://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__,
`native axes projections <https://matplotlib.org/stable/api/projections_api.html>`__,
`axis scales <https://matplotlib.org/stable/gallery/scales/scales.html>`__,
`colormaps <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`__,
`box styles <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyBboxPatch.html>`__,
`arrow styles <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyArrowPatch.html>`__,
and `arc styles <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.ConnectionStyle.html>`__
are referenced with "registered" string names,
as are `basemap projections <https://matplotlib.org/basemap/users/mapsetup.html>`__.
So, why not "register" everything else?

.. rubric:: Solution

In ProPlot, tick locators, tick formatters, axis scales, property cycles, colormaps,
normalizers, and `cartopy`_ projections are all "registered". This is accomplished
by defining "constructor functions" and passing various keyword arguments through
these functions.

The constructor functions also accept intuitive inputs alongside "registered"
names. For example, a scalar passed to `~proplot.constructor.Locator`
returns a `~matplotlib.ticker.MultipleLocator`, a
lists of strings passed to `~proplot.constructor.Formatter` returns a
`~matplotlib.ticker.FixedFormatter`, and `~proplot.constructor.Cycle`
and `~proplot.constructor.Colormap` accept colormap names, individual colors, and
lists of colors. Passing the relevant class instance to a constructor function
simply returns it, and all the registered classes are available in the top-level
namespace -- so class instances can be directly created with e.g.
``pplt.MultipleLocator(...)`` or ``pplt.LogNorm(...)`` rather than
relying on constructor functions.

For details, see the user guide sections on :ref:`Cartesian plots <ug_cartesian>`,
:ref:`color cycles <ug_cycles>` and :ref:`colormaps <ug_cmaps>`. The below table
lists the constructor functions and the keyword arguments that use them. Note that
`~matplotlib.axes.Axes.set_xscale` and `~matplotlib.axes.Axes.set_yscale` accept
instances of `~matplotlib.scale.ScaleBase` thanks to a patch applied by ProPlot.

================================  ============================================================  =============================================================  =================================================================================================================================================================================================
Function                          Return type                                                   Used by                                                        Keyword argument(s)
================================  ============================================================  =============================================================  =================================================================================================================================================================================================
`~proplot.constructor.Proj`       `~cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap`  `~proplot.ui.subplots`                                         ``proj=``
`~proplot.constructor.Locator`    `~matplotlib.ticker.Locator`                                  `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``locator=``, ``xlocator=``, ``ylocator=``, ``minorlocator=``, ``xminorlocator=``, ``yminorlocator=``, ``ticks=``, ``xticks=``, ``yticks=``, ``minorticks=``, ``xminorticks=``, ``yminorticks=``
`~proplot.constructor.Formatter`  `~matplotlib.ticker.Formatter`                                `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``formatter=``, ``xformatter=``, ``yformatter=``, ``ticklabels=``, ``xticklabels=``, ``yticklabels=``
`~proplot.constructor.Scale`      `~matplotlib.scale.ScaleBase`                                 `~proplot.axes.Axes.format`                                    ``xscale=``, ``yscale=``
`~proplot.constructor.Colormap`   `~matplotlib.colors.Colormap`                                 :ref:`2D plotting commands <ug_2dplots>`                        ``cmap=``
`~proplot.constructor.Norm`       `~matplotlib.colors.Normalize`                                :ref:`2D plotting commands <ug_2dplots>`                        ``norm=``
`~proplot.constructor.Cycle`      `~cycler.Cycler`                                              :ref:`1D plotting commands <ug_1dplots>`                        ``cycle=``
================================  ============================================================  =============================================================  =================================================================================================================================================================================================

.. _why_spacing:

Automatic dimensions and spacing
================================

.. rubric:: Limitation

Matplotlib plots tend to require lots of "tweaking" when you have more than one
subplot in the figure. This is partly because you must specify the physical dimensions
of the figure, despite the fact that...

#. The *subplot* aspect ratio is generally more relevant than the figure
   aspect ratio. An aspect ratio of ``1`` is desirable for most plots, and
   the aspect ratio must be held fixed for
   :ref:`geographic and polar <ug_proj>` projections and most
   `~matplotlib.axes.Axes.imshow` plots.
#. The physical width and height of the *subplot* controls the "evident"
   thickness of text, lines, and other content plotted inside the subplot.
   The effect of the figure size on this "evident" thickness depends on the
   number of subplot tiles in the figure.

Also, while matplotlib's `tight layout algorithm
<https://matplotlib.org/stable/tutorials/intermediate/tight_layout_guide.html>`__
can help you avoid tweaking the *spacing*, the algorithm cannot apply different
amounts of spacing between different subplot row and column boundaries.

.. rubric:: Solution

In ProPlot, you can specify the physical dimensions of a *reference subplot* by
passing `refwidth`, `refheight`, and/or `refaspect` to `~proplot.ui.figure` or
`~proplot.ui.subplots`. The dimensions of the figure are calculated automatically.
The default behavior is ``refaspect=1`` and ``refwidth=2`` (inches). If the
`aspect ratio mode
<https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axis_equal_demo.html>`__
for the reference subplot is set to ``'equal'``, as with :ref:`geographic <ug_geo>`,
:ref:`polar <ug_polar>`, `~matplotlib.axes.Axes.imshow`, and `~proplot.axes.Axes.heatmap`
plots, the :ref:`or data aspect ratio is used instead <ug_autosize>`.

You can also independently specify the width or height of the *figure* with the
`figwidth` and `figheight` parameters. If only one is specified, the other will be
calculated to preserve subplot aspect ratios. You can also select a `figwidth` and/or
`figheight` suitable for submission to :ref:`various publications <journal_table>`
using the `journal` parameter.

ProPlot also uses :ref:`its own tight layout algorithm <ug_tight>` to
automatically determine the `left`, `right`, `bottom`, `top`, `wspace`, and
`hspace` spacing parameters. This algorithm has the following advantages:

* Spacing between rows and columns is *variable* thanks to the
  `proplot.gridspec.GridSpec` subclass of `matplotlib.gridspec.GridSpec`.
  This is critical for putting :ref:`colorbars and legends <ug_cbars_legends>` or
  :ref:`axes panels <ug_insets_panels>` outside of subplots
  without "stealing space" from the parent subplot.
* The tight layout calculations are considerably simplified by permitting
  only *one* `proplot.gridspec.GridSpec` per figure. This restriction is
  possible by requiring successive `~proplot.figure.Figure.add_subplot`
  calls to imply the same geometry and include only subplot specs
  generated from the same `~proplot.gridspec.GridSpec`.

See the :ref:`user guide <ug_layout>` for details.

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


.. _why_redundant:

Working with multiple subplots
==============================

.. rubric:: Limitation

When working with multiple subplots in matplotlib, the path of least resistance
often leads to *redundant* figure elements. Namely...

* Repeated axis tick labels.
* Repeated axis labels.
* Repeated colorbars.
* Repeated legends.

These sorts of redundancies are very common even in publications, where they waste
valuable page space. It is also generally necessary to add "a-b-c" labels to
figures with multiple subplots before submitting them to publications, but
matplotlib has no built-in way of doing this.

.. rubric:: Solution

ProPlot makes it easier to work with multiple subplots and create clear,
concise figures.

* Axis tick labels and axis labels are automatically
  :ref:`shared and aligned between subplots <ug_share>` in the same
  `~proplot.gridspec.GridSpec` row or column. This is controlled by the `sharex`,
  `sharey`, `spanx`, `spany`, `alignx`, and `aligny` figure keywords.
* The figure `proplot.figure.Figure.colorbar` and `proplot.figure.Figure.legend`
  commands can easily draw colorbars and legends intended to reference more than
  one subplot in arbitrary contiguous rows and columns. See the
  :ref:`next section <why_colorbars_legends>` for details.
* The `~proplot.axes.Axes.panel_axes` (shorthand `~proplot.axes.Axes.panel`) commands
  can draw :ref:`thin panels <ug_panels>` along the edges of subplots. This
  can be useful for plotting 1D summary statistics alongside 2D plots.
* :ref:`A-b-c labels <ug_abc>` can be added to subplots simply using the :rcraw:`abc`
  setting -- for example, ``pplt.rc['abc'] = 'A.'`` or ``axs.format(abc='A.')``.
  This is possible because `~proplot.figure.Figure.add_subplot` assigns a unique
  `~proplot.axes.Axes.number` to every added subplot.
* The `proplot.figure.SubplotGrid.format` command can easily format multiple subplots
  at once or add colorbars, legends, panels, twin axes, or inset axes to multiple
  subplots at once. A `~proplot.figure.SubplotGrid` is returned by
  `proplot.figure.Figure.subplots`, and can be indexed like a list or like a 2D
  array (in which case the indices match the subplot grid extents).
  See the :ref:`user guide <ug_subplotgrid>` for details.


.. _why_colorbars_legends:

Simpler colorbars and legends
=============================

.. rubric:: Limitation

In matplotlib, it can be difficult to draw `~matplotlib.figure.Figure.legend`\ s
along the outside of subplots. Generally, you need to position the legend
manually and tweak the spacing to make room for the legend.

Also, `~matplotlib.figure.Figure.colorbar`\ s drawn along the outside of subplots
with e.g. ``fig.colorbar(..., ax=ax)`` need to "steal" space from the parent subplot.
This can cause asymmetry in figures with more than one subplot. It is also generally
difficult to draw "inset" colorbars in matplotlib and to generate outer colorbars
with consistent widths (i.e., not too "skinny" or "fat").

.. rubric:: Solution

ProPlot includes a simple framework for drawing colorbars and legends
that reference :ref:`individual subplots <ug_cbars_axes>` and
:ref:`multiple contiguous subplots <ug_cbars_figure>`.

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
rows and columns, "outer" colorbars and legends do not alter subplot
spacing or add whitespace. This is critical e.g. if you have a
colorbar between columns 1 and 2 but nothing between columns 2 and 3.
Also, `~proplot.figure.Figure` and `~proplot.axes.Axes` colorbar widths are
now specified in *physical* units rather than relative units, which makes
colorbar thickness independent of subplot size and easier to get just right.

There are also several useful new :ref:`colorbar <ug_cbars>` and
:ref:`legend <ug_legends>` features described in the user guide.


.. _why_plotting:

Improved plotting commands
==========================

.. rubric:: Limitation

A few common plotting tasks take a lot of work using matplotlib alone. The `seaborn`_,
`xarray`_, and `pandas`_ packages offer improvements, but it would be nice to
have this functionality built right into matplotlib's interface.

.. rubric:: Solution

ProPlot uses the `~proplot.axes.PlotAxes` subclass to add various `seaborn`_,
`xarray`_, and `pandas`_ features to existing matplotlib plotting commands
along with several additional features designed to make your life easier.

The following features are relevant for the 1D plotting commands like
`~proplot.axes.PlotAxes.line` (equivalent to `~proplot.axes.PlotAxes.plot`)
and `~proplot.axes.PlotAxes.scatter`:

* The `cycle` keyword is interpreted by the `~proplot.constructor.Cycle`
  :ref:`constructor function <why_constructor>` and applies
  :ref:`property cyclers <ug_apply_cycle>` on-the-fly. This permits succinct
  and flexible property cycler declaration.
* The `legend` and `colorbar` keywords draw :ref:`on-the-fly legends and colorbars
  <ug_cbars_axes>` using the result of the plotting command. Note that colorbars can
  be drawn from :ref:`lists of artists <ug_cbars>` (see `~proplot.axes.Axes.legend`).
* The default `ylim` (`xlim`) in the presence of a fixed `xlim` (`ylim`) is now
  adjusted to exclude out-of-bounds data. This can be useful when "zooming in" on
  a dependent variable axis but can be disabled by setting :rcraw:`axes.inbounds`
  to ``False`` or passing ``inbounds=False`` to plot commands.
* The `~proplot.axes.PlotAxes.bar` and `~proplot.axes.PlotAxes.barh` commands accept 2D
  arrays and can :ref:`stack or group <ug_bar>` successive columns. Likewise, the
  `~proplot.axes.PlotAxes.area` and `~proplot.axes.PlotAxes.areax` commands (shorthands
  for `~proplot.axes.PlotAxes.fill_between` and `~proplot.axes.PlotAxes.fill_betweenx`)
  accept 2D arrays and can :ref:`stack or overlay <ug_bar>` successive columns.
* The `~proplot.axes.PlotAxes.bar`, `~proplot.axes.PlotAxes.barh`,
  `~proplot.axes.PlotAxes.vlines`, `~proplot.axes.PlotAxes.hlines`,
  `~proplot.axes.PlotAxes.area`, and `~proplot.axes.PlotAxes.areax`
  commands accept a `negpos` keyword argument that assigns different
  colors to "negative" and "positive" regions.
* The `~proplot.axes.PlotAxes.linex` and `~proplot.axes.PlotAxes.scatterx` commands
  are just like `~proplot.axes.PlotAxes.line` and `~proplot.axes.PlotAxes.scatter`,
  but positional arguments are interpreted as *x* coordinates or (*y*, *x*) pairs.
  There are also the related commands `~proplot.axes.PlotAxes.stemx`,
  `~proplot.axes.PlotAxes.stepx`, `~proplot.axes.PlotAxes.boxh` (shorthand for
  `~proplot.axes.PlotAxes.boxploth`), and `~proplot.axes.PlotAxes.violinh` (shorthand
  for `~proplot.axes.PlotAxes.violinploth`).
* The `~proplot.axes.PlotAxes.line`, `~proplot.axes.PlotAxes.linex`,
  `~proplot.axes.PlotAxes.scatter`, `~proplot.axes.PlotAxes.scatterx`,
  `~proplot.axes.PlotAxes.bar`, and `~proplot.axes.PlotAxes.barh` commands can
  quickly draw vertical or horizontal :ref:`error bars or "shading" <ug_errorbars>`
  using a variety of keyword arguments. This is often more convenient than
  working directly with `~matplotlib.axes.Axes.errorbar`.
* The `~proplot.axes.PlotAxes.parametric` command draws clean-looking
  :ref:`parametric lines <ug_parametric>` by encoding the parametric
  coordinate using colormap colors rather than text annotations.

The following features are relevant for the 2D plotting commands like
`~proplot.axes.PlotAxes.pcolor` and `~proplot.axes.PlotAxes.contour`:

* The `cmap` and `norm` :ref:`keyword arguments <ug_apply_cmap>` are interpreted
  by the `~proplot.constructor.Colormap` and `~proplot.constructor.Norm`
  :ref:`constructor functions <why_constructor>`. This permits succinct
  and flexible colormap and normalizer application.
* The `colorbar` keyword draws on-the-fly :ref:`colorbars <ug_cbars_axes>`
  using the result of the plotting command. Note that "inset" colorbars can also
  be drawn, analogous to "inset" legends (see `~proplot.axes.Axes.colorbar`).
* The `~proplot.axes.PlotAxes.contour`, `~proplot.axes.PlotAxes.contourf`,
  `~proplot.axes.PlotAxes.pcolormesh`, and `~proplot.axes.PlotAxes.pcolor` commands
  all accept a `labels` keyword. This draws :ref:`contour and grid box labels
  <ug_labels>` on-the-fly. Labels are automatically colored black or white
  according to the luminance of the underlying grid box or filled contour.
* The default `vmin` and `vmax` used to normalize colormaps now excludes data
  outside the *x* and *y* axis bounds `xlim` and `ylim` if they were explicitly
  fixed. This can be disabled by setting :rcraw:`image.inbounds` to ``False``
  or by passing ``inbounds=False`` to plot commands.
* The `~proplot.colors.DiscreteNorm` normalizer is paired with most colormaps by
  default. It can easily divide colormaps into distinct levels, similar to contour
  plots. This can be disabled by setting :rcraw:`image.discrete` to ``False`` or
  by passing ``discrete=False`` to plot commands.
* The `~proplot.colors.DivergingNorm` normalizer is perfect for data with a
  :ref:`natural midpoint <ug_norm>` and offers both "fair" and "unfair" scaling.
  The `~proplot.colors.SegmentedNorm` normalizer can generate
  uneven color gradations useful for :ref:`unusual data distributions <ug_norm>`.
* The `~proplot.axes.PlotAxes.heatmap` command invokes
  `~proplot.axes.PlotAxes.pcolormesh` then applies an :ref:`equal axes apect ratio
  <https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axis_equal_demo.html>`,
  adds ticks to the center of each gridbox, and disables minor ticks and gridlines.
  This can be convenient for things like covariance matrices.
* Coordinate centers passed to commands like `~proplot.axes.PlotAxes.pcolor` are
  automatically translated to "edges", and coordinate edges passed to commands
  like `~proplot.axes.PlotAxes.contour` are automatically translated to "centers".
  In matplotlib, ``pcolor`` simply truncates the data when it receives centers.
* Commands like `~proplot.axes.PlotAxes.pcolor`, `~proplot.axes.PlotAxes.contourf`
  and `~proplot.axes.Axes.colorbar` automatically fix an irritating issue where
  saved vector graphics appear to have thin white lines between `filled contours
  <https://stackoverflow.com/q/8263769/4970632>`__, `grid boxes
  <https://stackoverflow.com/q/27092991/4970632>`__, and `colorbar segments
  <https://stackoverflow.com/q/15003353/4970632>`__. This can be disabled by
  passing ``fixedges=False`` to plot commands.

.. _why_cartopy_basemap:

Cartopy and basemap integration
===============================

.. rubric:: Limitation

There are two widely-used engines for working with geographic data in
matplotlib: `cartopy`_ and `basemap`_.  Using cartopy tends to be
verbose and involve boilerplate code, while using basemap requires plotting
with a separate `~mpl_toolkits.basemap.Basemap` object rather than the
`~matplotlib.axes.Axes`. They both require separate import statements and extra
lines of code to configure the projection.

Furthermore, when you use `cartopy`_ and `basemap`_ plotting
commands, "map projection" coordinates are the default coordinate system
rather than longitude-latitude coordinates. This choice is confusing for
many users, since the vast majority of geophysical data are stored with
longitude-latitude (i.e., "Plate Carrée") coordinates.

.. rubric:: Solution

ProPlot can succinctly create detailed geographic plots using either cartopy
or basemap as "backends". By default, cartopy is used, but basemap can be used
by passing ``basemap=True`` to axes-creation commands or by setting :rcraw:`basemap`
to ``True``. To create a geographic plot, simply pass the `PROJ <https://proj.org>`__
name to an axes-creation command, e.g. ``fig, ax = pplt.subplots(proj='pcarree')``
or ``fig.add_subplot(proj='pcarree')``. Alternatively, use the
`~proplot.constructor.Proj` constructor function to quickly generate
a `cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap` instance.

Requesting geographic projections creates a `proplot.axes.GeoAxes`
with unified support for `cartopy`_ and `basemap`_ features via the
`proplot.axes.GeoAxes.format` command. This lets you quickly modify geographic
plot features like latitude and longitude gridlines, gridline labels, continents,
coastlines, and political boundaries. The syntax is conveniently analogous to the
syntax used for `proplot.axes.CartesianAxes.format` and `proplot.axes.PolarAxes.format`.

The `~proplot.axes.GeoAxes` subclass also makes longitude-latitude coordinates
the "default" coordinate system by passing ``transform=ccrs.PlateCarree()``
or ``latlon=True`` to plotting commands (depending on whether cartopy or basemap
is the backend). And to enforce global coverage over the poles and across longitude
seams, you can pass ``globe=True`` to 2D plotting commands like
`~proplot.axes.PlotAxes.contour` and `~proplot.axes.PlotAxes.pcolormesh`.

See the :ref:`user guide <ug_proj>` for details.


.. _why_xarray_pandas:

Pandas and xarray integration
=============================

.. rubric:: Limitation

Scientific data is commonly stored in array-like containers
that include metadata -- namely, `xarray.DataArray`\ s, `pandas.DataFrame`\ s,
and `pandas.Series`. When matplotlib receives these objects, it simply ignores
the associated metadata. To create plots that are labeled with the metadata,
you must use the `xarray.DataArray.plot`, `pandas.DataFrame.plot`,
and `pandas.Series.plot` commands instead.

This approach is fine for quick plots, but not ideal for complex ones. It requires
learning a different syntax from matplotlib, and tends to encourage using the
`~matplotlib.pyplot` interface rather than the object-oriented interface. The
``plot`` commands also include features that would be useful additions to matplotlib
in their own right, without requiring special containers and a separate interface.

.. rubric:: Solution

ProPlot reproduces many of the `xarray.DataArray.plot`,
`pandas.DataFrame.plot`, and `pandas.Series.plot` features on the
`~proplot.axes.Axes` plotting commands themselves.  Passing a
`~xarray.DataArray`, `~pandas.DataFrame`, or `~pandas.Series` through any
plotting command updates the axis tick labels, axis labels, subplot title, and
colorbar and legend labels from the metadata. This feature can be disabled
by setting :rcraw:`autoformat` to ``False`` or passing ``autoformat=False``
to any plotting command.

ProPlot also supports `pint.Quantity` positional arguments by auto-calling
`~pint.UnitRegistry.setup_matplotlib` when a `pint.Quantity` is detected and by
extracting magnitudes from *z* coordinates (e.g., the data passed to ``contour``)
to avoid the stripped-units warning message. It also adds a unit string formatted
with :rcraw:`unitformat` as the default *x* and *y* axis label when :rcraw:`autoformat`
is enabled and supports `~xarray.DataArray` containers with `pint.Quantity` arrays.

Finally, as :ref:`described above <why_plotting>`, ProPlot implements features
that were originally only available from the `xarray.DataArray.plot`,
`pandas.DataFrame.plot`, and `pandas.Series.plot` commands -- like grouped
bar plots, layered area plots, and on-the-fly colorbars and legends --
directly within the `~proplot.axes.Axes` plotting commands.


.. _why_aesthetics:

Aesthetic colors and fonts
==========================

.. rubric:: Limitation

A common problem with scientific visualizations is the use of "misleading"
colormaps like ``'jet'``. These colormaps have jarring jumps in
`hue, saturation, and luminance <rainbow_>`_ that can trick the human eye into seeing
non-existing patterns. It is important to use "perceptually uniform" colormaps
instead. Matplotlib comes packaged with `a few of its own <matplotlib_>`_, plus
the `ColorBrewer <brewer_>`_ colormap series, but external projects
offer a larger variety of aesthetically pleasing "perceptually uniform" colormaps.

Matplotlib also "registers" the X11/CSS4 color names, but these are relatively
limited. The more intuitive and more numerous `XKCD color survey <xkcd_>`_ names
can be accessed with the ``'xkcd:'`` prefix, but this is cumbersome, and external
projects like `open color <opencolor_>`_ offer even more useful names.

Finally, matplotlib comes packaged with ``DejaVu Sans`` as the default font.
This font is open source and include glyphs for a huge variety of characters,
but unfortunately (in our opinion) it is not very aesthetically pleasing. It
can also be difficult to change the default matplotlib font.

.. rubric:: Solution

ProPlot adds new colormaps, colors, and fonts to help you make more
aesthetically pleasing figures.

* ProPlot adds colormaps from the `seaborn <seacolor_>`_, `cmocean <cmocean_>`_,
  `SciVisColor <sciviscolor_>`_, and `Scientific Colour Maps <fabio_>`_ projects.
  It also defines a few default :ref:`perceptually uniform colormaps <ug_perceptual>`
  and includes a `~proplot.colors.PerceptualColormap` class for generating
  new ones. A :ref:`table of colormap <ug_cmaps_included>` and
  :ref:`color cycles <ug_cycles_included>` can be shown using
  `~proplot.demos.show_cmaps` and `~proplot.demos.show_cycles`.
  Colormaps like ``'jet'`` can still be accessed, but this is discouraged.
* ProPlot adds colors from the `open color <opencolor_>`_ project and adds
  `XKCD color survey <xkcd_>`_ names without the ``'xkcd:'`` prefix after
  *filtering* them to exclude perceptually-similar colors and *normalizing* the
  naming pattern to make them more self-consistent. Old X11/CSS4 colors can still be
  accessed, but this is discouraged. A :ref:`table of color names <ug_colors_included>`
  can be shown using `~proplot.demos.show_colors`.
* ProPlot adds the entire `TeX Gyre <texgyre_>`_ font family to matplotlib. These
  are open-source fonts designed to resemble more popular, commonly-used fonts like
  Helvetica and Century. They are used as the new default serif, sans-serif, monospace,
  cursive, and "fantasy" fonts, and they are available on all workstations.
  A :ref:`table of font names <ug_fonts_included>` can be shown
  using `~proplot.demos.show_fonts`.

For details on adding new colormaps, colors, and fonts, see the
:ref:`.proplot folder <why_dotproplot>` section.

.. _why_colormaps_cycles:

Manipulating colormaps and cycles
=================================

.. rubric:: Limitation

In matplotlib, colormaps are implemented with the
`~matplotlib.colors.LinearSegmentedColormap` class (representing "smooth"
color gradations) and the `~matplotlib.colors.ListedColormap` class (representing
"categorical" color sets). They are generally cumbersome to modify or create from
scratch. Meanwhile, property cycles used for individual plot elements are implemented
with the `~cycler.Cycler` class. They are also cumbersome to modify and they cannot be
"registered" by name like colormaps.

The `seaborn`_ package introduces "color palettes" to make working with colormaps
and property cycles easier, but it would be nice to have similar features integrated
more closely with matplotlib.

..
   Colormap identification is also suboptimal, since the names are case-sensitive, and
   reversed versions of each colormap are not guaranteed to exist.

.. rubric:: Solution

In ProPlot, it is easy to manipulate colormaps and property cycles.

* All colormaps in ProPlot are replaced with the `~proplot.colors.ContinuousColormap`
  and `~proplot.colors.DiscreteColormap` subclasses of
  `~matplotlib.colors.LinearSegmentedColormap` and `~matplotlib.colors.ListedColormap`.
  These classes include several useful features leveraged by the
  :ref:`constructor functions <ug_constructor>`
  `~proplot.constructor.Colormap` and `~proplot.constructor.Cycle`.
* The `~proplot.constructor.Colormap` function can merge, truncate, and
  modify existing colormaps or generate brand new colormaps. It can also
  create new `~proplot.colors.PerceptualColormap`\ s -- a type of
  `proplot.colors.ContinuousColormap` with linear transitions in the
  :ref:`perceptually uniform-like <ug_perceptual>` hue, saturation,
  and luminance channels rather then the red, blue, and green channels.
* The `~proplot.constructor.Cycle` function can make property cycles from
  scratch or retrieve "registered" color cycles from their associated
  `~proplot.colors.DiscreteColormap` instances. It can also make property
  cycles by splitting up the colors from registered or on-the-fly
  `~proplot.colors.ContinuousColormap`\ s and `~proplot.colors.PerceptualColormap`\ s.

ProPlot also makes all colormap and color cycle names case-insensitive, and
colormaps are automatically reversed or cyclically shifted 180 degrees if you
append ``'_r'`` or ``'_s'`` to any colormap name. These features are powered by
`~proplot.colors.ColormapDatabase`, which replaces matplotlib's native
colormap database.

.. _why_norm:

Physical units engine
=====================

.. rubric:: Limitation

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
interpreting `figsize`, `figwidth`, `figheight`, `refwidth`, `refheight`,
`left`, `right`, `top`, `bottom`, `wspace`, `hspace`, and keyword arguments in a
few other places. Acceptable units include inches, centimeters, millimeters,
pixels, `points <https://en.wikipedia.org/wiki/Point_(typography)>`__, `picas
<https://en.wikipedia.org/wiki/Pica_(typography)>`__, and `em-heights
<https://en.wikipedia.org/wiki/Em_(typography)>`__ (a table of acceptable units
is found :ref:`here <units_table>`). Em-heights are particularly useful, as the
figure text can be a useful "ruler" when figuring out the amount of space you
need. The `~proplot.utils.units` function also translates rc settings assigned
to `~proplot.config.rc_matplotlib` and `~proplot.config.rc_proplot`, e.g.
:rcraw:`axes.labelpad`, :rcraw:`legend.handlelength`, and
:rcraw:`subplot.refwidth`. See the :ref:`user guide <ug_units>` for details.

.. _why_rc:

Flexible global settings
========================

.. rubric:: Limitation

In matplotlib, there are several `~matplotlib.rcParams` that would be
useful to set all at once, like spine and label colors. It might also
be useful to change these settings for individual subplots rather
than globally.

.. rubric:: Solution

In ProPlot, you can use the `~proplot.config.rc` object to change both native
matplotlib settings (found in `~proplot.config.rc_matplotlib`) and added ProPlot
settings (found in `~proplot.config.rc_proplot`). Assigned settings are always
validated, and special settings like ``meta.edgecolor``, ``meta.linewidth``, and
``font.smallsize`` can be used to update many settings all at once. Settings can
be changed with ``pplt.rc.key = value``, ``pplt.rc[key] = value``,
``pplt.rc.update(key=value)``, using `proplot.axes.Axes.format`, or using
`proplot.config.Configurator.context`. Settings that have changed during the
python session can be saved to a file with `proplot.config.Configurator.save`
(see `~proplot.config.Configurator.changed`), and settings can be loaded from
files with `proplot.config.Configurator.load`. See the
:ref:`user guide <ug_config>` for details.

.. _why_dotproplot:

Loading saved settings
======================

.. rubric:: Limitation

Matplotlib `~matplotlib.rcParams` can be changed persistently by placing
``matplotlibrc`` files in the same directory as your python script. But it
can be difficult to design and store your own colormaps and color cycles for
future use. It is also difficult to get matplotlib to use custom ``.ttf`` and
``.otf`` font files, which may be desirable when you are working on
Linux servers with limited font selections.

.. rubric:: Solution

ProPlot settings can be changed persistently by editing the default ``proplotrc``
file in the location given by `~proplot.config.Configurator.user_file` (this is
usually ``$HOME/.proplot/proplotrc``) or by adding ``proplotrc`` files to either
the current directory or any parent directory. Adding files to parent directories
can be useful when working in projects with lots of subfolders. See the
:ref:`user guide <ug_proplotrc>` for details.

ProPlot also automatically registers colormaps, color cycles, colors, and font
files stored in the ``cmaps``,  ``cycles``, ``colors``, and ``fonts`` folders in
the location given by `~proplot.config.Configurator.user_folder` (this is usually
``$HOME/.proplot``). You can save colormaps and color cycles to these
folders simply by passing ``save=True`` to `~proplot.constructor.Colormap` and
`~proplot.constructor.Cycle`. To explicitly register objects stored in these folders,
or to register arbitrary input arguments, you can use `~proplot.config.register_cmaps`,
`~proplot.config.register_cycles`, `~proplot.config.register_colors`, or
`~proplot.config.register_fonts`.
