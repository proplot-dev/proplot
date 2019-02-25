Figure and subplot panels
=========================

Outer Panels
------------

ProPlot figures may optionally have “outer” panels on the bottom, left,
or right of the figure, accessed with ``fig.bottompanel``,
``fig.leftpanel``, and ``fig.rightpanel``, respectively. They are instances of `~proplot.axes.PanelAxes`, and can be set up with `~proplot.subplots.subplots` in several ways:

-  ``bottompanel=True``: Allot space for a single panel spanning all
   columns of subplots.
-  ``bottompanels=True``: Allot space for ``n`` separate panels spanning
   the ``n`` columns of subplots.
-  ``bottompanels=[array]``: Allot space for an arbitrary number of
   panels spanning contiguous subplot columns – for example,
   ``bottompanels=[1,1,1,2,2,3]`` draws 3 panels, spanning 3 columns, 2
   columns, and 1 column respectively.

To access the nth panel, use ``fig.bottompanel[n]``. If you didn't request a panel in the call to `~proplot.subplots.subplots`, thi will raise a ``"Panel does not exist."`` attribute error.

If you call `~proplot.subplots.subplots` with ``bottomcolorbar[s]`` instead of ``bottompanel[s]``, ProPlot will modify the default axes widths to be *suitable for colorbars*. Otherwise, they default a bit thicker – intended for plotting stuff.

The `~proplot.axes.PanelAxes` differ from `~proplot.axes.XYAxes` in the `colorbar` and `legend` methods;

* ``fig.bottompanel.legend(handles, **kwargs)`` makes
  the panel axes invisible, then fills the space with a legend.
* ``fig.bottompanel.colorbar(mappable, **kwargs)`` turns the panel axes
  into a *colorbar*.

These are useful for global legends/colorbars intended to reference multiple subplots. Here's a simple example:

.. code-block:: python

   f, ax = plot.subplots(bottompanel=True, rightpanel=True)
   m = ax.contourf(np.random.rand(10,10))
   f.bototmpanel.colorbar(m, length=0.8)
   lines = ax.plot(np.random.rand(10,10)) # 10 lines, colored according to the active property cycler
   f.rightpanel.colorbar(lines, values=np.arange(10)) # see next section

Outer panel settings can be configured with a bunch of `subplots` keyword arguments; see the :ref:`documentation` for details.

Colorbar Enhancements
---------------------

Normally, to create a colorbar, you need a “mappable” instance – i.e.,
something with a `get_cmap` method, returned by `contourf`,
`pcolor`, etc.

With ProPlot, you have the following two additional options:

* ``colorbar(colors, values=values)`` creates a `~matplotlib.colors.ListedColormap` using a list of color strings or ``(R,G,B)`` tuples, then draws a colorbar with the values in the iterable `values` mapped to each color.
* ``colorbar(handles, values=values)`` infers colors from a list of “plot handles” -- i.e. anything with a `get_color` method, including handles returned by `~matplotlib.axes.Axes.plot`).

Two additional options make it easy to configure your colorbar geometry:

* ``length=fraction``, where ``0 <= fraction <= 1``, will make the colorbar span a *fraction* of the horizontal/vertical extent of the axes it is filling.
* ``extendlength=size`` controls the length of the “triangles” representing out-of-bounds colors (drawn when you use ``extend='min'``, ``extend='max'``, or ``extend='both'``). Since the “triangles” are now specified in physical units, they will always match other colorbars in the figure.

Inner Panels
------------

ProPlot also provides utilities for making “**inner panels**”. These may
be useful where you want a colorbar for every plot, a legend outside of
every axes, or want to show the x/y-direction statistics for some 2D
value plotted in your subplot (e.g. the x-direction mean, variance,
etc.).

The procedure for requesting inner panels is similar:

* ``innerpanels='r'`` draws panels on the right of each subplot, ``innerpanels='rt'`` draws panels on the top and the right.
* ``innerpanels={1:'r', (2,3):''}`` or ``innerpanels={range(5):'bt', 5:''}`` draws inner panels for *particular subplot numbers*.
* ``innerpanels_kw={'key':value}`` or, for example, ``innerpanels_kw={1:{'key':value1}, range(1,3):{'key':value2}}`` will format your inner panels.

As with the outer panels, you can also use ``innercolorbars='r'`` to
draw panels with default widths/spacing suitable for colorbars.

Inner panel settings can be configured with a bunch of `subplots` keyword arguments; see the :ref:`documentation` for details.

