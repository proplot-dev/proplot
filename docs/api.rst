.. _api:

=============
API reference
=============

The comprehensive API reference. All of the below objects are imported
into the top-level namespace. Use ``help(pplt.object)`` to read
the docs during a python session.

Please note that the documentation of plotting command "wrappers" from
ProPlot < 0.8 is now found in the individual `~proplot.axes.PlotAxes`
plotting commands. Using ``help(ax.command)`` during a python session shows both
the ProPlot documentation and the original matplotlib documentation.

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


Top-level functions
===================

.. automodule:: proplot.ui

.. automodsumm:: proplot.ui
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
