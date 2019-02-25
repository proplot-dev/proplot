Formatting your figures
=======================

The format command
------------------
The new `~proplot.axes.XYAxes.format` method, available on every axes returned from `~proplot.subplots.subplots`, is a versatile and powerful tool. But some might argue
that it just replicates features already available from `matplotlib` -- so, some motivation is in order.

To modify an axes property (e.g. an *x*-axis label) with the default API, you normally have to use a bunch of one-liner `~matplotlib.pyplot` commands, or method calls on axes and axis instances. This can get quite repetitive and quite verbose, resulting in lots of ugly, cumbersome boilerplate code.

Now, you can pass these settings to `~proplot.axes.XYAxes.format`. Instead of having to remember the name of the function, whether it's attached to the `~matplotlib.pyplot` module or an object instance, and the order and names of the arguments, you just have to remember one thing -- the name of the `~proplot.axes.XYAxes.format` argument.

The `~proplot.axes.XYAxes.format` method also abstracts away some inconsistencies and redundancies in the matplotlib API -- now, *There's Only One (obvious) Way To Do It*.

Example:

.. code-block:: python

   import proplot as plot
   f, ax = plot.subplots()
   ax.format(xlabel='time (seconds)', ylabel='temperature (K)', title='20th century sea-surface temperature')

Axis scales
-----------

This package adds several new axis scales, which can be invoked with
``[x|y]scale='name'`` in calls to ``format()``. They are described
below.

-  The **inverse** scale ``'inverse'``. Useful for, e.g., having
   wavenumber and wavelength on opposite sides of the same plot.
-  The **sine-weighted** and **Mercator** axis scales, ``'sine'`` and
   ``'mercator'``. The former creates an area-weighted latitude axis.
-  The configurable **cutoff** scale ``cutoff``. This allows arbitrary
   zoom-ins/zoom-outs over segments of an axis. Configure this scale
   using ``[x|y]scale=('cutoff', scale, lower, upper)`` where ``lower``
   and ``upper`` are the boundaries within which the axis scaling is
   multiplied by ``scale``. Use ``np.inf`` for a hard cutoff.
   Alternatively, use ``[x|y]scale=('cutoff', scale, position)`` to
   scale every coordinate after position ``position`` by ``scale`` (e.g.
   ``10`` would make the axis “faster”)

Tick locators
-------------

Control tick positions with ``[x|y]locator=arg`` or ``[x|y]ticks=arg``,
or supply it with a locator instance using ``plot.locator()``. The
locator can be specified in any of several ways:

-  Use e.g. ``xticks=5`` to tick every ``5`` data values.
-  Use ``xticks=[array]`` to tick specific locations.
-  Use ``xticks='string'`` to look up any of the ``matplotlib.ticker``
   locators, e.g. ``locator='month'`` or ``locator='log'``.

Finally, I recommend using ``plot.arange`` to generate lists of ticks –
it’s like ``np.arange``, but is **endpoint-inclusive**, which is
generally what you’ll want for this usage.

Tick formatters
---------------

Control tick label formatters with ``[x|y]formatter=arg`` or
``[x|y]ticklabels=arg``, or supply it with a formatter instance using
``plot.formatter()``. The new default ``CustomFormatter`` class for
ticklabels renders numbers into the style you’ll want 90% of the time.

I’ve also created several special formatter classes:

* Basic formatters ``'%.f'``, for classic `%-style
formatting <https://pyformat.info/>`__, or ``{x}`` for newer
``'string'.format(x=value)`` formatting.
* Coordinate formatters ``'lat'``, ``'deglat'``, ``'lon'``, ``'deglon'``, ``'deg'``: Use these
to format axis labels with cardinal direction indicators and/or degree
symbols, as denoted by the respective names.
* Fraction formatters ``'pi'``, ``'e'``, ``('symbol', scale)``: For tick labels represented as
fractions of some symbol.
* Explicit tick labels with e.g. ``xticklabels=['a', 'b', 'c']`` – adds specific strings to existing major ticks.

