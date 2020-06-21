Installation
============

ProPlot is published on `PyPi <https://pypi.org/project/proplot/>`__
and `conda-forge <https://conda-forge.org>`__. It can be installed
with ``pip`` or ``conda`` as follows:

.. code-block:: bash

   pip install proplot
   conda install -c conda-forge proplot

Likewise, an existing installation of ProPlot can be upgraded to the latest version with:

.. code-block:: bash

   pip install --upgrade proplot
   conda upgrade proplot


To install a development version of ProPlot, you can use
``pip install git+https://github.com/lukelbd/proplot.git``
or clone the repository and run ``pip install -e .`` inside
the ``proplot`` folder.

ProPlot's only hard dependency is `matplotlib <https://matplotlib.org/>`__.
The *soft* dependencies are `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__,
`basemap <https://matplotlib.org/basemap/index.html>`__,
`xarray <http://xarray.pydata.org>`__, and `pandas <https://pandas.pydata.org>`__.
See the documentation for details.
