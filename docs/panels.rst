Figure and subplot panels
=========================

Outer Panels
------------

ProPlot figures may optionally have “outer” panels on the bottom, left,
or right of the figure, accessed with ``fig.bottompanel``,
``fig.leftpanel``, and ``fig.rightpanel``, respectively. They can be set
up by providing any of several arguments to ``plot.subplots()`` (we use
``bottom`` as an example):

-  ``bottompanel=True``: Allot space for a single panel spanning all
   columns of subplots.
-  ``bottompanels=True``: Allot space for ``n`` separate panels spanning
   the ``n`` columns of subplots.
-  ``bottompanels=[array]``: Allot space for an arbitrary number of
   panels spanning contiguous subplot columns – for example,
   ``bottompanels=[1,1,1,2,2,3]`` draws 3 panels, spanning 3 columns, 2
   columns, and 1 column respectively.

To access the ``n``\ th panel, use ``fig.bottompanel[n]``. If you didn’t
request a panel in the call to ``plot.subplots()``, then try to access
it, you will get a ``"Panel does not exist."`` attribute error.

If you call ``plot.subplots()`` with ``bottomcolorbar[s]`` instead of
``bottompanel[s]`` (or ``leftcolorbar[s]`` instead of ``leftpanel[s]``,
``rightcolorbar[s]`` instead of ``rightpanel[s]``), ProPlot will modify
the default axes widths to be *suitable for colorbars* (otherwise, they
default a bit thicker – intended for plotting stuff).

Access axes plotting methods like normal, with
``fig.bottompanel.method`` or ``fig.bottompanel[n].method``, with two
caveats: \* Use ``fig.bottompanel.legend(handles, **kwargs)`` to make
that axes invisible, then fill the space with a legend. This is useful
for creating **global legends intended to reference multiple subplots**.
\* Use ``fig.bottompanel.colorbar(mappable, **kwargs)`` to turn that
axes into a *colorbar*.

Outer panel geometry can be configured with the following
``plot.subplots()`` keyword arguments. Note once again that, while
*numeric* sizing arguments are assumed to be in inches, you can specify
sizes with **arbitrary units** using e.g. ``width='12cm'``,
``wspace='5mm'``, ``hspace='3em'`` (3 em squares), etc. \* ``lwidth``,
``rwidth``, ``bwidth``: Widths of the left, right, and bottom panels,
respectively. \* ``lspace``, ``rspace``, ``bspace``: Empty space on the
*outside* (i.e. toward the figure edge) of left, right, and bottom
panels, respectively.

Example:

.. code:: python

   f, ax = plot.subplots(bottompanel=True, rightpanel=True)
   m = ax.contourf(np.random.rand(10,10))
   f.bototmpanel.colorbar(m, length=0.8)
   lines = ax.plot(np.random.rand(10,10)) # 10 lines, colored according to the active property cycler
   f.rightpanel.colorbar(lines, values=np.arange(10)) # see next section

Colorbar Enhancements
---------------------

Normally, to create a colorbar, you need a “mappable” instance – i.e.,
something with a ``get_cmap()`` method, returned by ``contourf``,
``pcolor``, etc.

With ProPlot, you have the following two additional options: \*
``colorbar([color-strings-or-tuples], values=values``: Create a
``ListedColormap`` using a list of color strings or ``(R,G,B)`` tuples,
then draw a colorbar where the values in the iterable ``values`` are
mapped to each color. \* ``colorbar([plot-handles], values=[values])``:
As above, but this time the colors are inferred from a list of “plot
handles” (anything with a ``get_color()`` method, e.g. the handles
returned by ``ax.plot()``).

Two additional options make it easy to configure your colorbar geometry:
\* Use ``length=fraction`` for ``0 <= fraction <= 1`` to make the
colorbar span a *fraction* of the horizontal/vertical extent of the axes
it is filling. \* Use ``extendlength=size`` in *inches* to control the
length of the “triangles” representing out-of-bounds colors (drawn when
you use ``extend='min'``, ``extend='max'``, or ``extend='both'``). The
“triangles” are now always specified in **physical units**, so they will
look consistent with other colorbars in the figure.

Inner Panels
------------

ProPlot also provides utilities for making “**inner panels**”. These may
be useful where you want a colorbar for every plot, a legend outside of
every axes, or want to show the x/y-direction statistics for some 2D
value plotted in your subplot (e.g. the x-direction mean, variance,
etc.).

The procedure for requesting inner panels is similar: \* Use
``innerpanels='r'`` to draw inner subplot panels on the right. Use
``innerpanels='rt'`` to draw one panel on the right, another on the top.
\* Use dictionary arguments, for example
``innerpanels={1:'r', (2,3):''}`` or
``innerpanels={range(5):'bt', 5:''}``, to draw inner panels for only
**particular subplot numbers**. \* Format your inner panels with
``innerpanels_kw={'key':value}`` or, for example,
``innerpanels_kw={1:{'key':value1}, range(1,3):{'key':value2}}`` to
format the panels differently.

The inner panel dimensions can be specified with the following keyword
arguments, using ``innerpanels_kw``: \* ``wwidth``, ``hwidth``: Width of
vertical (left/right) and horizontal (top/bottom) panels, respectively.
\* ``wspace``, ``hspace``: Empty space on the *inside* (toward the main
subplot) of left, right, and bottom panels, respectively. Use e.g.
``lspace=0`` to join the panel with the main subplot.

As with the outer panels, you can also use ``innercolorbars='r'`` to
draw panels with default widths/spacing suitable for colorbars.
