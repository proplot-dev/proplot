==============
Zen of ProPlot
==============

ProPlot is meant to improve upon the parts of matplotlib that
tend to be laborious or cumbersome
when trying to create
beautiful graphics. This page
enumerates these limitations and
describes how ProPlot addresses them.

Efficiency
==========
* Power users often need to change lots of plot settings all at once.
* In matplotlib, this requires a bunch of one-liner setters and getters, like `~matplotlib.axes.Axes.set_title`. It is often unclear whether settings can be changed from the `~matplotlib.axes.Axes` instance or the ``ax.xaxis`` and ``ax.yaxis`` axis instances, or whether an ordinary setter exists or something more complex, e.g. `~matplotlib.axes.Axes.tick_params`, must be used. In general, this workflow is cumbersome and results in boilerplate plotting code.
* In ProPlot, the `~proplot.axes.Axes.format` command is introduced for changing settings in bulk.

For example, it is trivial to see that

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
* Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, `~cartopy.crs.AlbersEqualArea`, and `~cartopy.crs.NorthPolarStereo`. With the default APIs, these classes must be invoked directly.
* With ProPlot, these classes are all *registered*, just like `axis scales <https://matplotlib.org/3.1.0/gallery/scales/scales.html>`__ and `colormaps <https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html>`__. This is done by creating *constructor functions* and passing various keyword argument through the functions, allowing users to reference "registered" names instead of the actual classes. The below table lists the functions and the keyword arguments that use them.

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

Handy axes containers
=====================
.. The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.figure.Figure.add_subplot` commands.

* In matplotlib, `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself. This variable output is cumbersome.
* In ProPlot, `~proplot.subplots.subplots` returns an `axes_grid` of axes that unifies the behavior of these three possible return values. `axes_grid` is a `list` subclass that behaves like a scalar when it contains just one element, supports row-major or column-major 1D indexing (e.g. ``axs[0]``), and permits 2D indexing (e.g. ``axs[1,0]``) no matter the geometry. Since `~proplot.subplots.subplots` can generate figures with arbitrarily complex subplot geometry, this 2D indexing is useful only when the arrangement happens to be a clean 2D matrix. Further, thanks to the `~axes_grid.__getattr__` override, the `~proplot.subplots.axes_grid` allows you to call arbitrary methods on arbitrary axes all at once, e.g. ``axs.format(tickminor=False)``.

Fluid figure dimensions
=======================
* In matplotlib, you have to specify the physical dimensions of the figure. However, the dimensions of the *individual subplots* are often more important, since (1) they determine subplot aspect ratios and (2) this controls the *evident* thickness of e.g. a 12 point font or 2 point line relative to the subplot box.
* In ProPlot, you can specify the physical size of *subplots* instead of the figure by passing ``axwidth`` or ``axheight`` to `~proplot.figure.Figure`, and the default behavior is ``axwidth=2``. Figure dimensions are then allowed to vary to accommodate the subplot geometry and the "tight layout" spacing adjustments.

It is likely that matplotlib does not implement "flexible dimensions" because several matplotlib backends require figure size to be fixed. When `~matplotlib.figure.Figure.draw` changes the figure size, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

The right layout every time
===========================
* In matplotlib, the tight layout algorithm is cumbersome.
* In ProPlot, the tight layout algorithm is more accurate because `~matplotlib.gridspec.GridSpec`\ s permit variable spacing between rows and columns, and their spacing parameters are specified in physical units instead of figure-relative units.

.. The `~matplotlib.gridspec.GridSpec` class is useful for creating figures with complex subplot geometry.
.. Users want to control axes positions with gridspecs.
.. * Matplotlib permits arbitrarily many `~matplotlib.gridspec.GridSpec`\ s per figure. This greatly complicates the tight layout algorithm for little evident gain.

To simplify the tight layout algorithm, ProPlot permits only *one* `~matplotlib.gridspec.GridSpec` per figure. When a `~matplotlib.gridspec.SubplotSpec` is passed to `~proplot.subplots.Figure.add_subplot`, the figure is locked to the associated `~matplotlib.gridspec.GridSpec`. When an integer or tuple is passed to `~proplot.subplots.Figure.add_subplot`, the geometry implied by subsequent integer or tuple calls must *divide* or *multiply* the initial geometry -- for example, two square subplots above a longer rectangle subplot can be drawn by passing the integers ``221``, ``222``, and ``212`` to `~proplot.subplots.Figure.add_subplot`. This introduces a marginal limitation (see discussion in #50) but *considerably* simplifies the tight layout algorithm. 

Arbitrary units
===============
.. * Configuring spaces and dimensions in matplotlib often requires physical units.

* Matplotlib uses "inches" for figure dimensions and figure-relative or axes-relative units almost everywhere else. However, "inches" are foreign to the world outside of the U.S., and *relative* units encourage "tinkering" with meaningless numbers that change the subjective appearance when the figure dimensions change, since *text* and *plotted content* are specified in the physical units "points".
* ProPlot permits arbitrary physical units for almost all sizing arguments, e.g. ``left='0.5cm'``. This prevents "tinkering" and encourages users to be aware of the physical dimensions describing their figure. It also introduces font-relative units, e.g. ``left='1em'``, which are generally more intuitive than figure-relative units. Finally, you can return to axes-relative and figure-relative units with e.g. ``left='0.1fig'`` or ``left='0.1ax'``.

Working with colormaps
======================
* In matplotlib, colormaps are hard to modify and hard to create. Also, colormap names are case-sensitive, and reversed ``'_r'`` versions are not guaranteed to exist.
* In ProPlot, it is easy to generate, combine, and modify colormaps using the `~proplot.styletools.Colormap` constructor function. Also, all colormap names are case-insensitive, and all colormaps are reversible by appending ``'_r'`` to the name. "Colormaps" and "color cycles" are also *fluid*, e.g. you can use a colormap as the color cycler for line plots. This is ProPlot's answer to seaborn's "palettes".

.. -- matplotlib's "colormaps" and "property cyclers" are sufficient.

More accurate colorbars
=======================
* In matplotlib, when ``extend!='both'``, excess colors are trimmed from either side of the map. However, most of the time, users want the full color range of the colormap. Matplotlib usually divides colormaps up into levels by sampling the colormap with a low resolution lookup table, which makes it hard to get boundaries just right.
* In ProPlot, colorbars always traverse the full range of colormap colors, no matter the ``extend`` setting. ProPlot samples the colormap with a high resolution lookup table and uses `~proplot.styletools.BinNorm` to restrict the possible colormap indices. It makes more sense to handle *all* aspects of the value --> color conversion process, including *discretization*, to the normalizer.

New font families
=================
* In matplotlib, the default font is DejaVu Sans. This font is not very aesthetically pleasing.
* In ProPlot, the default font is Helvetica. This is a tried and tested, aesthetically pleasing sans serif font.

Cartopy and basemap defaults
============================
* In basemap and cartopy, the default coordinate system is always *map projection coordinates*. This is inexplicable and confusing for new users, whose data is almost always in longitude-latitude or "Plate Carrée" coordinates.
* In ProPlot, `latlon=True` is the basemap default, and `transform=ccrs.PlateCarree()` is the cartopy default.
