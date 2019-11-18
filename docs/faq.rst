==========================
Frequently asked questions
==========================

Why do inline figures look different / why is the font so small?
================================================================

This question needs some background: These days, most publications prefer plots saved as **vector graphics** (i.e. a PDF, EPS, or SVG file) rather than **bitmap graphics** (i.e. PNG and JPG). Vector graphics are infinitely scalable, and the file sizes are often smaller than bitmap graphics. We recommend using vector graphics even when your plots are not destined for publication.

While bitmap graphics record data in pixels, vector graphics record data in physical units, e.g. inches and points (equal to 1/72 inches). When you save a vector graphic, the content sizes should be appropriate for embedding the plot in a document. For example, if an academic journal recommends 8-point font for plots, you should use 8-point font in your plotting code.

This is where matplotlib comes in:

* The default matplotlib backends produce really small, low-resolution, artifact-plagued jpeg figures. To make the backend-generated images legible, matplotlib uses a fairly large default figure width of 6.5 inches (usually only suitable for multi-panel plots) and a fairly large default font size of 10 points (most journals recommend 5-9 points).
* Unfortunately, these sizes often mean your figure has to be downscaled. And when a vector graphic is downscaled, the sizes you used used for fonts and lines in your plotting code are not the sizes that appear in the document. This defeats the purpose of saving the figure as a vector graphic!

To remedy this, ProPlot uses a **higher resolution** default inline backend, a **smaller** default font size, and makes the default figure size dependent on the **number of subplots** in the figure. ProPlot also adds a ``journal`` keyword argument to `~proplot.subplots.Figure`, which lets you apply figure dimensions from a particular academic journal standard.

With ProPlot, you should never have to shrink your figure when embedding it in a document. The physical sizes recorded in the file will be true.

.. users to enlarge their figure dimensions and font sizes so that content inside of the inline figure is visible -- but when saving the figures for publication, it generally has to be shrunk back down!


How does ProPlot differ from other matplotlib wrappers?
=======================================================

There is already a great matplotlib wrapper called `seaborn <https://seaborn.pydata.org/>`__, and `pandas <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html>`__ and `xarray <http://xarray.pydata.org/en/stable/plotting.html>`__ both offer convenient matplotlib plotting commands. What makes this project different?

While parts of ProPlot were inspired by these projects (in particular, ``rctools.py`` and ``colortools.py`` are modeled after seaborn), the goals are largely different:

* ProPlot, seaborn, pandas, and xarray all offer tools for generating rigid, simple, nice-looking plots from data stored in `~pandas.DataFrame`\ s and `~xarray.DataArray`\ s (ProPlot tries to apply labels from these objects, just like pandas and xarray).
* Unlike seaborn, pandas, and xarray, ProPlot also permits a *high an extremely high level of customization*, implements these features *arbitrary complex grids of subplots*, and includes powerful tools for working with *geographic* data.
* ProPlot reproduces many of the seaborn tools for working with color *and* introduces some really awesome new tools, like `~proplot.styletools.Colormap` and `~proplot.styletools.PerceptuallyUniformColormap`.
* ProPlot is built *right into the matplotlib API*, thanks to special subclasses of the `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes, while seaborn, pandas, and xarray are meant to be used separately from the matplotlib API. This massively reduces the learning curve for new ProPlot users.

In summary, this project is intended to unify the convenience of seaborn, pandas, and xarray plotting with the power and customizability of the underlying matplotlib API.

Why not contribute to matplotlib directly?
==========================================

Since ProPlot is built right into the matplotlib API, you might be wondering why we didn't contribute to the matplotlib project directly. The main answer is *speed* and *autonomy*, but there are a few practical limitations:

* Certain features directly conflict with matplotlib. For example, ProPlot's "smart tight layout" conflicts with matplotlib's `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ by permitting *fluid figure dimensions*. Also, the new `~proplot.subplots.GridSpec` class permits *variable spacing* between rows and columns, and specifies spacing in *physical units* instead of figure-relative and axes-relative units.
* Other features may be too redundant. For example, `~proplot.axes.Axes.format` is convenient, but the same tasks can be accomplished with existing axes and axis "setter" methods. Also, some of the functionality of `~proplot.subplots.subplots` can be replicated with `axes_grid1 <https://matplotlib.org/mpl_toolkits/axes_grid1/index.html>`__. Following `TOOWTDI <https://wiki.python.org/moin/TOOWTDI>`__ philosophy, these features should probably not be integrated.

Nevertheless, if any core matplotlib think that some of ProPlot's features should be added to matplotlib, please contact core developer `Luke Davis <mailto:lukelbd@gmail.com>`__ and let him know!
