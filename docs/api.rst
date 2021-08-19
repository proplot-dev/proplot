.. _api:

=============
API reference
=============

Comprehensive documentation of ProPlot functions and classes. All of these
objects are imported into the top-level namespace, so you can read the
documentation within python sessions using ``help(pplt.function_or_class)``.
Please note that the "wrapper" function documentation from proplot < 0.8
is now located on the individual plotting commands under
`proplot.axes.PlotAxes`. When calling ``help(axes.command)`` on
plotting commands during a python session, both the ProPlot
documentation and the original matplotlib documentation are shown.

Top-level functions
===================

.. automodule:: proplot.ui

.. automodsumm:: proplot.ui
   :toctree: api


Figure class
============

.. automodule:: proplot.figure

.. automodsumm:: proplot.figure
   :toctree: api
   :skip: SubplotsContainer


Gridspec class
==============

.. automodule:: proplot.gridspec

.. automodsumm:: proplot.gridspec
   :toctree: api


Axes classes
============

.. automodule:: proplot.axes

.. automodsumm:: proplot.axes
   :toctree: api


Configuration tools
===================

.. automodule:: proplot.config

.. automodsumm:: proplot.config
   :toctree: api
   :skip: inline_backend_fmt, RcConfigurator


Constructor functions
=====================

.. automodule:: proplot.constructor

.. automodsumm:: proplot.constructor
   :toctree: api
   :skip: Colors


Locators and formatters
=======================

.. automodule:: proplot.ticker

.. automodsumm:: proplot.ticker
   :toctree: api


Axis scale classes
==================

.. automodule:: proplot.scale

.. automodsumm:: proplot.scale
   :toctree: api


Colormaps and normalizers
=========================

.. automodule:: proplot.colors

.. automodsumm:: proplot.colors
   :toctree: api
   :skip: ListedColormap, LinearSegmentedColormap, PerceptuallyUniformColormap, LinearSegmentedNorm


Projection classes
==================

.. automodule:: proplot.crs

.. automodsumm:: proplot.crs
   :toctree: api


Miscellaneous tools
===================

.. automodule:: proplot.utils

.. automodsumm:: proplot.utils
   :toctree: api
   :skip: shade, saturate


Demo functions
==============

.. automodule:: proplot.demos

.. automodsumm:: proplot.demos
   :toctree: api
