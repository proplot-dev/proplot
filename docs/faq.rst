==========================
Frequently asked questions
==========================

What makes this project different?
==================================

There is already a great matplotlib wrapper called
`seaborn <https://seaborn.pydata.org/>`__. Also, `pandas
<https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.pplt.html>`__
and `xarray <http://xarray.pydata.org/en/stable/plotting.html>`__
both offer convenient matplotlib plotting commands.
How does ProPlot compare against these tools?

* ProPlot, seaborn, pandas, and xarray all offer tools for generating rigid, simple,
  nice-looking plots from data stored in `~pandas.DataFrame`\ s and
  `~xarray.DataArray`\ s (ProPlot tries to apply labels from these objects, just like
  pandas and xarray).
* ProPlot is integrated with *cartopy* and *basemap*. You will find plotting geophysical
  data in ProPlot to be much more concise than working with cartopy and basemap
  directly.
* ProPlot *expands upon* the seaborn tools for working with color and global settings.
  For example, see `~proplot.constructor.Colormap`,
  `~proplot.colors.PerceptualColormap`, and `~proplot.config.Configurator`.
* ProPlot *expands upon* matplotlib by fixing various quirks, developing a more
  advanced automatic layout algorithm, simplifying the process of drawing outer
  colorbars and legends, and much more.
* ProPlot is *built right into the matplotlib API*, thanks to special subclasses of the
  `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes, while seaborn,
  pandas, and xarray are meant to be used separately from the matplotlib API.

In a nutshell, ProPlot is intended to *unify the convenience of seaborn, pandas, and
xarray plotting with the power and customizability of the underlying matplotlib API*.

..
  So while ProPlot includes similar tools, the scope and goals are largely different.
  Indeed, parts of ProPlot were inspired by these projects -- in particular,
  ``setup.py`` and ``colortools.py`` are modeled after seaborn. However the goals and
  scope of ProPlot are largely different:

Why didn't you add to matplotlib directly?
==========================================

Since ProPlot is built right into the matplotlib API, you might be wondering why we
didn't contribute to the matplotlib project directly.

* Certain features directly conflict with matplotlib. For example, ProPlot's tight
  layout algorithm conflicts with matplotlib's `tight layout
  <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ by
  permitting *fluid figure dimensions*, and the new `~proplot.gridspec.GridSpec` class
  permits *variable spacing* between rows and columns and uses *physical units* rather
  than figure-relative and axes-relative units.
* Certain features are arguably too redundant. For example, `~proplot.axes.Axes.format`
  is convenient, but the same tasks can be accomplished with existing axes and axis
  "setter" methods. Also, some of the functionality of `~proplot.ui.subplots` can be
  replicated with `axes_grid1
  <https://matplotlib.org/mpl_toolkits/axes_grid1/index.html>`__. Following `TOOWTDI
  <https://wiki.python.org/moin/TOOWTDI>`__ philosophy, these features should probably
  not be integrated.

..
   * ProPlot design choices are made with the academic scientist working with ipython
     notebooks in mind, while matplotlib has a much more diverse base of hundreds of
     thousands of users. Matplotlib developers have to focus on support and API
     consistency, while ProPlot can make more dramatic improvements.

..
   Nevertheless, if any core matplotlib developers think that some of ProPlot's features
   should be added to matplotlib, please contact
   `Luke Davis <https://github.com/lukelbd>`__ and let him know!

Why do my inline figures look different?
========================================

These days, most publications prefer plots saved as
`vector graphics <https://en.wikipedia.org/wiki/Vector_graphics>`__ [1]_
rather than `raster graphics <https://en.wikipedia.org/wiki/Raster_graphics>`__ [2]_.
When you save vector graphics, the content sizes should be appropriate for embedding the
plot in a document (for example, if an academic journal recommends 8-point font for
plots, you should use 8-point font in your plotting code).

Most of the default matplotlib backends make low-quality, artifact-plagued jpegs. To
keep them legible, matplotlib uses a fairly large default figure width of 6.5 inches
(usually only suitable for multi-panel plots) and a slightly large default font size of
10 points (where most journals recommend 5-9 points). This means your figures have to be
downscaled so the sizes used in your plotting code are *not* the sizes that appear in
the document.

ProPlot helps you get your figure sizes *correct* for embedding them as vector graphics
inside publications.  It uses a slightly smaller default font size, calculates the
default figure size from the number of subplot rows and columns, and adds the `journal`
keyword argument to `~proplot.figure.Figure` which can be used to employ figure
dimensions from a particular journal standard.  To keep the inline figures legible,
ProPlot also employs a *higher quality* default inline backend.

.. [1] `Vector graphics <https://en.wikipedia.org/wiki/Vector_graphics>`__ use physical
   units (e.g. inches, `points <https://en.wikipedia.org/wiki/Point_(typography)>`__),
   are infinitely scalable, and often have much smaller file sizes than bitmap graphics.
   You should consider using them even when your plots are not destined for publication.
   PDF, SVG, and EPS are the most common formats.

.. [2] `Raster graphics <https://en.wikipedia.org/wiki/Raster_graphics>`__ use pixels
   and are *not* infinitely scalable. They tend to be faster to display and easier
   to view, but they are discouraged by most academic publishers. PNG and JPG are the
   most common formats.

..
   users to enlarge their figure dimensions and font sizes so that content inside of the
   inline figure is visible -- but when saving the figures for publication, it generally
   has to be shrunk back down!
