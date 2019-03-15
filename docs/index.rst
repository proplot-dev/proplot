.. proplot documentation master file, created by
   sphinx-quickstart on Wed Feb 20 01:31:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======
ProPlot
=======

An `object-oriented <https://matplotlib.org/api/api_overview.html>`__ `matplotlib <https://matplotlib.org/>`__ wrapper
that can help you make beautiful, publication-quality graphics.

Installation
============

This package is a work-in-progress. Currently there is no formal release
on PyPi (but it's coming soon!). For the time being, you may install directly from Github using:

.. code-block:: bash

   pip install git+https://github.com/lukelbd/proplot.git#egg=proplot

The dependencies are `matplotlib <https://matplotlib.org/>`_ and `numpy <http://www.numpy.org/>`_.  The optional geographic mapping features require `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ and/or `basemap <https://matplotlib.org/basemap/index.html>`_.

I recommend importing with

.. code-block:: python

   import proplot as plot

to differentiate ProPlot from the usual "``plt``" abbreviation used for the `~matplotlib.pyplot` module. If you are in an ipython notebook, I also recommend adding

.. code-block:: python

   plot.nbsetup()

The `~proplot.notebook.nbsetup` command enables two useful iPython extensions and greatly increases the resolution of your inline figures -- ProPlot applies relatively small font size and line width defaults, so that the physical sizes match what might appear in print. See the `~proplot.rcmod` documentation for more on how to change these defaults.

Overview
========

On import, a bunch of new colormaps and string color names are registered.
If this is all you want, and you don't care about other features, simply
import ProPlot at the top of your script. See `~proplot.colortools` for details, and the :ref:`Table of colormaps`, :ref:`Table of color cycles`, and :ref:`Table of colors`.

The remaining features mostly derive from the `~proplot.subplots.subplots` command, inspired
by the pyplot `~matplotlib.pyplot.subplots` command. This generates a scaffolding
of axes with new, helpful methods. It also has a bunch of other useful features, like "panels"
and built-in geographic projections.

The next most important utility is the ``format`` method on `~proplot.axes.BaseAxes`, which calls so-called ``smart_update`` methods on the `~proplot.axes.BaseAxes`, `~proplot.axes.XYAxes`, `~proplot.axes.CartopyAxes`, and `~proplot.axes.BasemapAxes` axes types. The latter three types can all be returned by `~proplot.subplots.subplots`, depending on the arguments you used. Use `~proplot.axes.BaseAxes.format` to fine-tune your axis properties, titles, labels, limits, and much more. See :ref:`Table of projections` for the projections available to cartopy and basemap.

To get started, check out the :ref:`General introduction`.
Hopefully, you will find this API to be **less verbose** and **more powerful** than the builtin `pyplot and object-oriented <https://matplotlib.org/api/api_overview.html>`__ matplotlib APIs.

.. This is just so top-level headers in the showcase.rst
   file appear as *subsections* in the documentation to
   a parent 'showcase' section.

.. Below replaces traditional autodoc invocation with
   .. automodule:: mods.set.tests
      :members:
      :show-inheritance:
   Use 'sphinx-apidoc ../proplot -o .' to auto-generate modules.rst

.. toctree::
   :maxdepth: 4
   :caption: Contents

   tutorial
   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
