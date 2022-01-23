.. _api:

=============
API reference
=============

The comprehensive API reference. All of the below objects are imported
into the top-level namespace. Use ``help(pplt.object)`` to read
the docs during a python session.

Please note that proplot removes the associated documentation when functionality
is deprecated (see :ref:`What's New <whats_new>`). However, proplot adheres to
`semantic versioning <https://semver.org>`__, which means old code that uses deprecated
functionality will 1) still work and 2) issue warnings rather than errors until the
first major release (i.e. version 1.0.0).

.. important::

   The color transformation functions like `to_rgba` and `scale_luminance` from
   proplot < 0.10.0 can now be found as methods on the new `~proplot.colors.Color`
   class. Please see `~proplot.colors.Color` for the relevant documentation. Note that
   old code that uses commands like ``pplt.to_rgba()`` and ``pplt.scale_luminance()``
   will still work -- but it will result in a deprecation warning.

.. important::

   The documentation for "wrapper" functions like `standardize_1d` and `cmap_changer`
   from proplot < 0.8.0 like can now be found under individual `~proplot.axes.PlotAxes`
   methods like `~proplot.axes.PlotAxes.pcolor`. Please see `~proplot.axes.PlotAxes` for
   the relevant documentation. Note that calling ``help(ax.method)`` in a python session
   will show both the proplot documentation and the original matplotlib documentation.

Figure class
============

.. automodule:: proplot.figure

.. automodsumm:: proplot.figure
   :toctree: api


Grid classes
============

.. automodule:: proplot.gridspec

.. automodsumm:: proplot.gridspec
   :toctree: api
   :skip: SubplotsContainer


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


Demo functions
==============

.. automodule:: proplot.demos

.. automodsumm:: proplot.demos
   :toctree: api


Miscellaneous functions
=======================

.. automodule:: proplot.utils

.. automodsumm:: proplot.utils
   :toctree: api
   :skip: shade, saturate
