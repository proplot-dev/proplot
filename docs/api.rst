.. _api:

=============
API reference
=============

The comprehensive API reference. All of the below objects are imported
into the top-level namespace. Use ``help(pplt.object)`` to read
the docs during a python session.

Please note that proplot removes the associated documentation when functionality
is deprecated (see :ref:`What's New <whats_new>`). However, proplot adheres to
`semantic versioning <https://semver.org>`__, which means old code that uses
deprecated functionality will still work and issue warnings rather than errors
until the first major release (version 1.0.0).

.. important::

   The color transformation functions like `to_rgba` and `scale_luminance` from
   proplot < 0.10.0 can now be found as methods on the new `~proplot.colors.Color`
   class. Note that old code that uses commands like ``pplt.to_rgba()`` and
   ``pplt.scale_luminance()`` will still work (but result in a deprecation warning).

.. important::

   The documentation for "wrapper" functions like `standardize_1d` and `cmap_changer`
   from proplot < 0.8.0 can now be found under individual `~proplot.axes.PlotAxes`
   methods like `~proplot.axes.PlotAxes.plot` and `~proplot.axes.PlotAxes.pcolor`. Note
   that calling ``help(ax.method)`` in a python session will show both the proplot
   documentation and the original matplotlib documentation.

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

.. automodule:: proplot.proj

.. automodsumm:: proplot.proj
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
   :skip: shade, saturate get_colors, set_hue, set_saturation, set_luminance, set_alpha, shift_hue, scale_saturation, scale_luminance, to_hex, to_rgb, to_xyz, to_rgba, to_xyza
