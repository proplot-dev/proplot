.. proplot documentation master file, created by
   sphinx-quickstart on Wed Feb 20 01:31:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======
ProPlot
=======

A library providing helpful and versatile plotting utilities
for making beautiful, publication-quality graphics with ease.

Installation
============

This package is a work-in-progress. Currently there is no formal release
on PyPi. However, feel free to install directly from Github using:

.. code-block:: bash

   pip install git+https://github.com/lukelbd/proplot.git#egg=proplot

Dependencies are `matplotlib` and `numpy`. The geographic mapping
mapping features require `cartopy` and/or `basemap`. I recommend importing
with

.. code-block:: python

   import proplot as plot

to differentiate ProPlot from the usual `plt` abbreviation for `~matplotlib.pyplot`.
If you are in an ipython notebook, I also recommend adding

.. code-block:: python

   plot.nbsetup()

This enables a couple useful extensions and greatly increases the resolution of your inline figures.
See `~proplot.notebook.nbsetup` for details.

Overview
========

On import, a bunch of new colormaps and string color names are registered.
If this is all you want, and you don't care about other features, simply
import ProPlot at the top of your script. See `~proplot.colortools` for details,
and the visualizations in :ref:`Smooth colormaps`, :ref:`Discrete colormaps`,
:ref:`New color names`.

The remaining features mostly derive from the `~proplot.subplots.subplots` command, inspired
by the pyplot `~matplotlib.pyplot.subplots` command. This generates a scaffolding
of axes with new, helpful methods. It also has a bunch of other useful features, like "panels"
and built-in geographic projections.

The next most important utilities are the `format` methods on the `~proplot.axes.BaseAxes`,
`~proplot.axes.XYAxes`, `~proplot.axes.CartopyAxes`, and `~proplot.axes.BasemapAxes` axes
returned by `~proplot.subplots.subplots`. Use `format` to fine-tune
your axis properties, titles, labels, limits, and much more. This API is **much less verbose**
and **much more powerful** than the builtin
`pyplot and object-oriented <https://matplotlib.org/api/api_overview.html>`__
matplotlib APIs.

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

   header
   quickstart
   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
