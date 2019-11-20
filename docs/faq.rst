==========================
Frequently asked questions
==========================

How does ProPlot differ from other matplotlib wrappers?
=======================================================

There is already a great matplotlib wrapper called `seaborn <https://seaborn.pydata.org/>`__. Furthermore, `pandas <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html>`__ and `xarray <http://xarray.pydata.org/en/stable/plotting.html>`__ both offer convenient matplotlib plotting commands. How does ProPlot compare against these tools?

* ProPlot, seaborn, pandas, and xarray all offer tools for generating rigid, simple, nice-looking plots from data stored in `~pandas.DataFrame`\ s and `~xarray.DataArray`\ s (ProPlot tries to apply labels from these objects, just like pandas and xarray).
* Unlike seaborn, pandas, and xarray, (1) ProPlot includes new tools that permit an extremely high level of customization, (2) its features work for arbitrarily complex subplot grids instead of just single-subplot figures, and (3) it includes powerful tools for working with *geographic* data.
* ProPlot reproduces many of the seaborn tools for working with color and global settings, but also introduces some really awesome *new* tools, like `~proplot.styletools.Colormap`, `~proplot.styletools.PerceptuallyUniformColormap`, and `~proplot.rctools.rc_configurator`.
* ProPlot is built *right into the matplotlib API*, thanks to special subclasses of the `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes, while seaborn, pandas, and xarray are meant to be used separately from the matplotlib API. This massively reduces the learning curve for new ProPlot users.

In a nutshell, ProPlot is intended to *unify the convenience of seaborn, pandas, and xarray plotting with the power and customizability of the underlying matplotlib API*.

.. So while ProPlot includes similar tools, the scope and goals are largely different.
.. Indeed, parts of ProPlot were inspired by these projects -- in particular, ``rctools.py`` and ``colortools.py`` are modeled after seaborn. However the goals and scope of ProPlot are largely different:

Why not contribute to matplotlib directly?
==========================================

Since ProPlot is built right into the matplotlib API, you might be wondering why we didn't contribute to the matplotlib project directly. The main answer is *speed* and *autonomy*, but there are a few practical limitations:

* Certain features directly conflict with matplotlib. For example, ProPlot's "smart tight layout" conflicts with matplotlib's `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ by permitting *fluid figure dimensions*. Also, the new `~proplot.subplots.FlexibleGridSpec` class permits *variable spacing* between rows and columns, and specifies spacing in *physical units* instead of figure-relative and axes-relative units.
* Other features may be too redundant. For example, `~proplot.axes.Axes.format` is convenient, but the same tasks can be accomplished with existing axes and axis "setter" methods. Also, some of the functionality of `~proplot.subplots.subplots` can be replicated with `axes_grid1 <https://matplotlib.org/mpl_toolkits/axes_grid1/index.html>`__. Following `TOOWTDI <https://wiki.python.org/moin/TOOWTDI>`__ philosophy, these features should probably not be integrated.

Nevertheless, if any core matplotlib think that some of ProPlot's features should be added to matplotlib, please contact core developer `Luke Davis <https://github.com/lukelbd>`__ and let him know!

Why do inline figures look different / why is the font so small?
================================================================
These days, most publications prefer plots saved as `vector graphics <https://en.wikipedia.org/wiki/Vector_graphics>`__ [1]_ (e.g. PDF, EPS, SVG) rather than `raster graphics <https://en.wikipedia.org/wiki/Raster_graphics>`__ (e.g. PNG, JPG). The latter records data in pixels, while the former records data in physical units, e.g. inches and `points <https://en.wikipedia.org/wiki/Point_(typography)>`__. When you save a vector graphic, the content sizes should be appropriate for embedding the plot in a document. For example, if an academic journal recommends 8-point font for plots, you should use 8-point font in your plotting code.

This is where matplotlib comes in:

* The default matplotlib backends produce really small, low-resolution, artifact-plagued jpeg figures. To make the backend-generated images legible, matplotlib uses a fairly large default figure width of 6.5 inches (usually only suitable for multi-panel plots) and a fairly large default font size of 10 points (most journals recommend 5-9 points).
* Unfortunately, these sizes often mean your figure has to be downscaled. And when a vector graphic is downscaled, the sizes you used for fonts and lines in your plotting code are not the sizes that appear in the document. This kind of defeats the purpose of vector graphics.

To remedy this, ProPlot uses a **higher resolution** default inline backend, a **smaller** default font size, and makes the default figure size dependent on the **number of subplots** in the figure. ProPlot also adds a ``journal`` keyword argument to `~proplot.subplots.Figure`, which lets you apply figure dimensions from a particular academic journal standard.

With ProPlot, you should never have to shrink your figure when embedding it in a document. The physical sizes recorded in the file will be true.

.. [1] `Vector graphics <https://en.wikipedia.org/wiki/Vector_graphics>`__ are infinitely scalable, and the file sizes are often smaller than bitmap graphics. We recommend using vector graphics even when your plots are not destined for publication.

.. users to enlarge their figure dimensions and font sizes so that content inside of the inline figure is visible -- but when saving the figures for publication, it generally has to be shrunk back down!


