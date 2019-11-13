============
Why ProPlot?
============

ProPlot is meant to improve upon the parts of matplotlib that
tend to be laborious or cumbersome
when trying to create
beautiful graphics. This page
enumerates these limitations and
describes how ProPlot addresses them.

Efficiency
==========
* Power users generally want to change a ton of plot settings rather than relying on defaults.
* In matplotlib, this requires a bunch of one-liner setters and getters, like `~matplotlib.axes.Axes.set_title`. This is cumbersome and generally results in copy-paste plotting code.
* In ProPlot, we encourage changing settings in bulk, and introduce the `~proplot.axes.Axes.format` command for doing so. For example, it is trivial to see that

  ```python
  import proplot as plot
  f, ax = plot.subplots(ncols=2)
  ax.format(xticks=20, xtickminor=False, xlabel='x axis', ylabel='y axis')
  ```

  is much more succinct than

  ```python
  import matplotlib.pyplot as plt
  fig, ax = plt.subplots()
  ax.set_xticks(np.arange(0,100,10))
  ax.minorticks_off()
  ax.set_xlabel('x axis')
  ax.set_ylabel('y axis')
  ```

Constructor functions
=====================
* Matplotlib and cartopy introduce a bunch of classes with verbose names like
`~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, `~cartopy.crs.AlbersEqualArea`, and `~cartopy.crs.NorthPolarStereo`. With the default APIs, these classes must be invoked directly.
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

Flexible dimensions
===================
* In matplotlib, you have to specify the physical of the figure. However, the physical dimensions of the *individual subplots* are often more important, since (1) they determine subplot aspect ratios and (2) this controls the relative appearance of e.g. a 12 point text or 2 point thick line looks.
* In ProPlot, you can specify the physical size of *subplots* instead of the figure. Figure dimensions are then allowed to vary to accommodate the subplot geometry and the "tight layout" spacing adjustments. This is not possible in matplotlib -- likely, the reason this wasn't implemented is several matplotlib backends require figure size to be fixed. When `~matplotlib.figure.Figure.draw` changes the figure size, it is "too late" for the backend and weird bugs can arise. ProPlot fixes these issues with the static notebook inline backend and the Qt popup backend (thanks to a detail of its implementation; this might break someday). However, these issues are fixable for the "notebook" and "macosx" backends, and possibly other popup backends.

Axes grids
==========
.. The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.figure.Figure.add_subplot` commands.
* Matplotlib's `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself. This variable output is cumbersome.
* ProPlot's `~proplot.subplots.subplots` returns an `axes_grid` of axes, meant to unify these three possible return values.  `axes_grid` is a `list` subclass supporting 1D indexing (e.g. ``axs[0]``), but permits 2D indexing (e.g. ``axs[1,0]``) *just in case* the user *happened* to draw a clean 2D matrix of subplots. The `~axes_grid.__getattr__` override also means it no longer matters whether you are calling a method on an axes or a singleton `axes_grid` of axes. Finally, `axes_grid` lets `subplots` support complex arrangements of subplots -- just use 1D indexing when they don't look like a 2D matrix.

Flexible units
==============
.. * Configuring spaces and dimensions in matplotlib often requires physical units.
* Matplotlib just uses "inches" for figure dimensions and relative units almost everywhere else. However, "inches" are foreign to the world outside of the U.S., and relative units encourage "tinkering" with meaningless numbers that change the subjective appearance when figure or axes dimensions change. This is because *text* and *plotted content* are specified in the physical units "points".
* ProPlot permits arbitrary physical units for almost all sizing arguments, e.g. ``left=0.5cm`` or ``left=2em``. This encourages users to always be aware of their physical units, and font-size relative coordinates are usually more intuitive than figure-relative coordinates. You can also return to axes-relative and figure-relative coordinates with e.g. ``left=0.1fig`` or ``left=0.1ax``.

Inline backend behavior
=======================
* The default inline ipython notebook backend produces really small, low-resolution, artifact-plagued jpeg figures. This encourages users to enlarge their figure dimensions and font sizes so that content inside of the inline figure is visible -- but when saving the figures for publication, it generally has to be shrunk back down!
* ProPlot uses a higher resolution default inline backend and smaller, more realistic default font and figure sizes. This means users won't usually have to downsize their figures when saving them as a PDF, SVG, or EPS vector graphics and embedding them in a PDF for publication -- the "points" and "inches" recorded in the vector graphics are *correct*.

One gridspec
============
.. The `~matplotlib.gridspec.GridSpec` class is useful for creating figures with complex subplot geometry.
.. Users want to control axes positions with gridspecs.

* Matplotlib permits arbitrarily many gridspecs per figure. This massively complicates tight layout algorithms and other operations for little evident gain.
* ProPlot permits only *one* gridspec per figure. When `subplots` is used, this is trivial to enforce. When `~proplot.subplots.Figure.add_subplot` is used, the figure geometry is "locked" after the first call -- although `~proplot.subplots.Figure.add_subplot` calls that divide into the existing geometry are also acceptable (for example, two square subplots above a longer rectangle subplots with the integers ``221``, ``222``, and ``212``). This choice is not a major imposition on the user (see discussion in #50), and *considerably* simplifies gridspec adjustments, e.g. the "tight layout" adjustments. Also, the ProPlot tight layout algorithm is more accurate because GridSpecs can now have variable spacing between rows and columns.

Nitty gritty colorbars
======================
* In matplotlib, when ``extend!='both'``, excess colors are trimmed from either side of the map. However, most of the time, time users want the full color range of the colormap. Matplotlib usually divides colormaps up into levels by sampling the colormap with a low resolution lookup table, which makes it hard to get boundaries just right.
* In ProPlot, colorbars always traverse the full range of colormap colors, no matter the ``extend`` setting. ProPlot samples the colormap with a high resolution lookup table and uses `~proplot.styletools.BinNorm` to restrict the possible colormap indices. IMO it makes more sense to handle *all* aspects of the value --> color conversion process, including *discretization*, to the normalizer.

Colormaps names and colormap conversion
=======================================
* In matplotlib, 
The color usage tools are incredibly useful. Should mention how colormap names are case-insensitive, all colormaps and cycles are reversible by appending `_r` to the name. Also, "cycles" and "colormaps" are *fluid*, e.g. you can use a colormap as the color cycler for line plots. This is ProPlot's answer to the seaboarn idea of "palettes", which I think are unnecessarily convoluted -- matplotlib's "colormaps" and "property cyclers" are sufficient.


Cartopy and basemap improvements
================================
* Better projection handling and more sensible defaults. Part of this falls into the "use constructor functions everywhere and register verbose class objects under succinct string names" philosophy.
* Another is that basemap `latlon=True` should **really** be the default (99% of the time people are not working with data stored in map projection coordinates!). Same goes with cartopy `transform=ccrs.PlateCarree()` -- I think cartopy made the map projection transform the default because they wanted to copy the behavior of basemap.
