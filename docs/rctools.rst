Global settings control
=======================

.. automodsumm:: proplot.rctools
   :toctree: api
   :skip: units, get_ipython

.. Hacky bullshit to prevent header appearing in TOC
.. raw:: html

   <h1>Summary</h1>

.. automodule:: proplot.rctools


.proplotrc file
---------------

You can edit the default global settings by placing a file in your home directory called ``.proplotrc``.
The ``.proplotrc`` file containing the ProPlot defaults is shown below. It roughly matches the syntax used in ``.matplotlibrc``
files, although we strictly adhere to `YAML <https://en.wikipedia.org/wiki/YAML>`__.

.. include:: ../proplot/.proplotrc
   :literal:
