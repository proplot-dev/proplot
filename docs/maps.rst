Geographic projections
======================
ProPlot also lets you set up axes with geographic projections using either of 2 packages: `~mpl_toolkits.basemap.Basemap` or `~cartopy.crs.Projection`. Projections are configured with the `proj` and `proj_kw` arguments for the `~proplot.subplots.subplots` command.

Control the map projection type with `proj='proj'` or e.g. `proj={1:'proj1', (2,3):'proj2', 4:'name3'}`. In the latter case, the integers and integer tuples correspond to **the subplot number**.

In the same way, you can pass keyword arguments (e.g. `lon_0`) to the cartopy `~cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap` classes using `proj_kw={'name':value}` or e.g. `proj_kw={1:'proj1', (2,3):'proj2'}`. You can also choose between cartopy and basemap using `basemap=False` or e.g. `basemap={1:True, 2:False}`.

Example:

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots(ncols=3, proj='hammer', basemap=False)

This creates 3 side-by-side `Hammer projections <https://en.wikipedia.org/wiki/Hammer_projection>`_ using Cartopy.

Note that `~mpl_toolkits.basemap` is no **`longer under active development <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`_** -- cartopy is the intended replacement, as it is integrated more intelligently with the matplotlib API.

However, for the time being, `~mpl_toolkits.basemap` *retains several advantages* over cartopy. Namely, `more tools for labeling meridians/parallels <https://github.com/SciTools/cartopy/issues/881>`_ and more available projections -- see the `basemap list <https://matplotlib.org/basemap/users/mapsetup.html>`_ vs. the `cartopy list <https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html>`_. I therefore decided to support both.

Cartopy axes
------------
When you specify the `proj` keyword arg with ``basemap=False``, a `~proplot.axes.CartopyAxes` instance (subclassed from the cartopy `~cartopy.geoaxes.Geoaxes` class) is created. As shown above, you can now declare the projection by **string name**, instead of having to reference cartopy `~cartopy.crs.Projection` classes directly.

In cartopy, you usually need to supply ``transform=cartopy.crs.PlateCarree()`` to the plotting method (see `this example <https://scitools.org.uk/cartopy/docs/v0.5/matplotlib/introductory_examples/03.contours.html>`_). With ProPlot, **this is done by default**.

Other aspects of cartopy axes can be controlled with `~proplot.axes.CartopyAxes.format`.

Basemap axes
------------
When you specify the `proj` keyword arg with ``basemap=True``, a `~proplot.axes.BasemapAxes` instance is created. This class allows you to access basemap plotting utilities **directly on the axes as a method**, instead of having to call the method from the `~mpl_toolkits.basemap.Basemap` instance.

To fix issues with the "seam" on the edge of the map, I've overridden several plotting methods on `~proplot.axes.BasemapAxes` -- data will be automatically circularly rolled until the left-hand-side comes after the map seam, then interpolated to the seam longitudes.

Other aspects of basemap axes can be controlled with `~proplot.axes.BaseAxes.format`.

