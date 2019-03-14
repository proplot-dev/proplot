Formatting
==========

The format command
------------------
The new `~proplot.axes.XYAxes.format` method, available on every axes returned from `~proplot.subplots.subplots`, is a versatile and powerful tool. But some might argue
that it just replicates features already available from matplotlib -- so, some motivation is in order.

To modify an axes property (e.g. an *x*-axis label) with the default API, you normally have to use a bunch of one-liner `~matplotlib.pyplot` commands, or method calls on axes and axis instances. This can get quite repetitive and quite verbose, resulting in lots of ugly, cumbersome boilerplate code.

Now, you can pass these settings to `~proplot.axes.XYAxes.format`. Instead of having to remember the name of the function, whether it's attached to the `~matplotlib.pyplot` module or an object instance, and the order and names of the arguments, you just have to remember one thing -- the name of the `~proplot.axes.XYAxes.format` argument.

The `~proplot.axes.XYAxes.format` method also abstracts away some inconsistencies and redundancies in the matplotlib API -- now, *There's Only One (obvious) Way To Do It*.

Example:

.. code-block:: python

   import proplot as plot
   f, ax = plot.subplots()
   ax.format(xlabel='time (seconds)', ylabel='temperature (K)', title='20th century sea-surface temperature')

Scale settings
--------------

Control axis scales with the `xscale` and `yscale` keyword args.
There are also some new scales available, described below:

-  The **inverse** scale ``'inverse'``. Useful for, e.g., having
   wavenumber and wavelength on opposite sides of the same plot.
-  The **sine-weighted** and **Mercator** axis scales, ``'sine'`` and
   ``'mercator'``. The former creates an area-weighted latitude axis.
-  The **cutoff** scale, allowing arbitrary
   zoom-ins and zoom-outs over segments of an axis. This scale is invoked
   with the tuple ``('cutoff', scale, lower, upper)`` where ``lower``
   and ``upper`` are the boundaries within which the axis scaling is
   multiplied by ``scale``. Use ``np.inf`` for a hard cutoff.
   Alternatively, use ``('cutoff', scale, position)`` to
   scale every coordinate after position ``position`` by ``scale``.

Tick settings
-------------

Control tick positions with the `xlocator` and `ylocator` keyword args (or their aliases, `xticks` and `yticks`). These accept a number of possible values:

*  A number (e.g. ``xticks=N``) ticks every N data values.
*  A string will look up any of the `matplotlib.ticker`
   locators by key name, e.g. ``'log'``.
*  A list of numbers will tick those specific locations.

I recommend using `plot.arange` to generate lists of ticks –
it’s like `numpy.arange`, but is **endpoint-inclusive**, which more often than
not is what you'll want in this context.

Control tick label formats with the `xformatter` and `yformatter` keyword args (or their aliases, `xticklabels` and `yticklabels`). These accept a number of possible values:

* ``'%.f'`` for classic `%-style formatting <https://pyformat.info/>`_, or ``{}`` for newer `'string'.format(value)` formatting.
* ``'lat'``, ``'deglat'``, ``'lon'``, ``'deglon'``, and ``'deg'``
  format axis labels with cardinal direction indicators and/or degree
  symbols, as denoted by the respective names.
* ``'pi'``, ``'e'``, ``('symbol', scale)`` will format tick labels represented as
  fractions of some symbol (the first 2 are :math:`\pi` and Euler's constant, provided for convenience).
* A list of strings (e.g. ``xticklabels=['a', 'b', 'c']``) will simply label existing ticks with that list.


The rc object
-------------
A special object named `~proplot.rcmod.rc` (belonging to a class called
`~proplot.rcmod.rc_configurator`) is created whenever you import ProPlot. This object
gives you advanced control over the look of your plots. **Use**
`~proplot.rcmod.rc` **as your one-stop shop for changing global settings**.

The `~proplot.rcmod.rc` object can be used to change built-in
`~matplotlib.rcParams` settings, a few custom :ref:`rcParams_new` settings,
and some magic :ref:`rcGlobals` settings that apply to groups of other
settings and keep them synced – e.g., tick, spine, and tick label
colors. The global settings are tabulated in the `~proplot.rcmod` documentation.

To modify any :ref:`rcGlobals` or `~matplotlib.rcParams` settings, you have four options:

1. Change the default settings for good by creating a `.proplotrc` file in your home folder. For more information, see the `~proplot.rcmod` documentation.
2. Change one global setting using ``plot.rc.name = value`` or ``plot.rc['name'] = value``.
   Note that, for settings with ‘dots’ in their name, you will
   have to use ``plot.rc['category.name'] = value``
3. Update several global settings at once using
   ``plot.rc.update({'name1':value1, 'name2':value2})`` or
   ``plot.rc.update(name1=value1, name2=value2)``, just like you would
   update a dictionary.
4. Change local settings using
   ``ax.format(rc_kw={'name1':value1, 'name2':value2})`` or
   ``ax.format(name1=value1, name2=value2)``. In this case, *the rc settings will only be applied to that specific axes*. This can be convenient for (e.g.) drawing focus to a particular subplot by changing
   its color. If the "rc" setting you want to change has a dot in its name, simply omit the dot -- for example, the custom ProPlot setting ``title.pos`` can be changed with ``ax.format(titlepos='ci')``.

To access a single setting, use ``rc.name`` or ``rc['name']``. To
access a group of setting by category name, use e.g. ``rc.axes``
and a dictionary of settings will be returned.

To reset everything to the default state, use `~proplot.rcmod.rc_configurator.reset`. By default, settings are reset every time a figure is drawn -- that is, when a figure is rendered by the matplotlib backend or saved to file.

