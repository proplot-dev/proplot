Formatting your figure
======================

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

Control axis scales with the `xscale` and `yscale` keyword args.
There are also some new scales available, described below:

-  The **inverse** scale ``'inverse'``. Useful for, e.g., having
   wavenumber and wavelength on opposite sides of the same plot.
-  The **sine-weighted** and **Mercator** axis scales, ``'sine'`` and
   ``'mercator'``. The former creates an area-weighted latitude axis.
-  The **cutoff** scale, allowing arbitrary
   zoom-ins/zoom-outs over segments of an axis. This scale is invoked
   with the tuple ``('cutoff', scale, lower, upper)`` where ``lower``
   and ``upper`` are the boundaries within which the axis scaling is
   multiplied by ``scale``. Use ``np.inf`` for a hard cutoff.
   Alternatively, use ``('cutoff', scale, position)`` to
   scale every coordinate after position ``position`` by ``scale`` (e.g.
   ``10`` would make the axis “faster”).

Tick settings
-------------

Control tick positions with the `xlocator` and `ylocator` keyword args (or their aliases, `xticks` and `yticks`). These accept a number of possible values:

*  A number (e.g. ``xticks=N``) ticks every ``5`` data values.
*  A string will look up any of the `matplotlib.ticker`
   locators by key name, e.g. ``'log'``.
*  A list of numbers will tick those specific locations.

I recommend using `plot.arange` to generate lists of ticks –
it’s like `np.arange`, but is **endpoint-inclusive**, which more often than
not is what you'll want in this context.

Control tick label formats with the `xformatter` and `yformatter` keyword args (or their aliases, `xticklabels` and `yticklabels`). These accept a number of possible values:

* ``'%.f'`` for classic `%-style formatting <https://pyformat.info/>`_, or ``{}`` for newer ``'string'.format(value)`` formatting.
* ``'lat'``, ``'deglat'``, ``'lon'``, ``'deglon'``, and ``'deg'``
  format axis labels with cardinal direction indicators and/or degree
  symbols, as denoted by the respective names.
* ``'pi'``, ``'e'``, ``('symbol', scale)`` will format tick labels represented as
  fractions of some symbol (the first 2 are :math:`\pi` and Euler's constant, provided for convenience).
* A list of strings (e.g. ``xticklabels=['a', 'b', 'c']``) will simply label existing ticks with that list.


The rc object
-------------
A special object named ``rc`` (belonging to a class called
`rc_configurator`) is created whenever you import ProPlot. This object
gives you advanced control over the look of your plots. **Use**
`plot.rc` **as your one-stop shop for changing global settings**.

``rc`` object can be used to change built-in
``matplotlib.rcParams`` settings, a few custom “``rcSpecial``” settings,
and some magic “``rcGlobal``” settings that apply to groups of other
settings and keep them synced – e.g., tick, spine, and tick label
colors.

The following is an overview of the available ``rcGlobal`` settings.

* ``color``, ``xcolor``, ``ycolor``: The color of axis spines, tick
  marks, tick labels, and labels. Use the ``x`` and ``y`` versions to just
  change settings for one axis.
* ``facecolor``, ``facehatch``: The
  background color and optional hatching pattern (e.g. ``'xxx'``; see
  `this demo <https://matplotlib.org/gallery/shapes_and_collections/hatch_demo.html>`__)
  for an axes. The latter is useful where you wish to highlight invalid
  (transparent) “NaN” data in a ``pcolormesh`` or ``contourf`` plot.
* ``small``, ``large``: Font size for legend text, tick labels, axis
  labels, and manually placed text created by ``ax.text`` (``small``), and
  for titles, super titles, and a-b-c labels (``large``).
* ``fontname``: The font name used for all text in the figure. ProPlot comes packaged
  with a bunch of desirable fonts, and **changes the default from DejaVu
  Sans (or Bitstream Vera) to Helvetica Neue**. When you first import
  ProPlot, run ``plot.install_fonts()`` and restart your ipython session
  to install them.
* ``linewidth``, ``minorwidth``: thickness of axes spines and major tick
  lines (``linewidth``), and minor tick lines (``minorwidth``, where minor tick line thickness = ``linewidth * minorwidth``).
* ``gridwidth``, ``gridratio``: thickness of gridlines (``gridwidth``), and
  minor gridlines (``gridratio``, where minor gridline thickness
  = ``gridwidth * gridratio``).
* ``gridalpha``, ``gridcolor``, ``gridstyle``: the transparency, color, and line style
  for your major and minor gridlines.
* ``ticklen``, ``tickratio``: length of major ticks (``ticklen``), and
  minor ticks (``tickratio``, where minor tick lengths = ``ticklen * tickratio``).
* ``tickdir``: tick direction, one of ``out``, ``in``, or ``inout``.
* ``abcweight``, ``titleweight``, ``suptitleweight``: the font weight (one of
  ``ultralight``, ``light``, ``normal``, ``medium``, ``demi``, ``bold``,
  ``very bold``, or ``black``) for your title text. Note that many fonts
  only have ``normal`` or ``bold`` available; if you request another
  weight, the “closest” availble weight will be selected.

To modify any ``rcGlobals`` or ``rcParams`` settings, you have three options:

1. Change one global setting using ``plot.rc.name = value`` or ``plot.rc['name'] = value``.
   Note that, for ``rcParams`` settings with ‘dots’ in their name, you will
   have to use ``plot.rc['category.name'] = value``.
2. Update several global settings at once using
   ``plot.rc.update({'name1':value1, 'name2':value2})`` or
   ``plot.rc.update(name1=value1, name2=value2)``, just like you would
   update a dictionary.
3. Change local settings using
   ``ax.format(rc_kw={'name1':value1, 'name2':value2})`` or
   ``ax.format(name1=value1, name2=value2)``. Note that, for this last
   option, **the rc settings will only be applied locally** (i.e. to the
   axes on which ``format()`` is being invoked). This can be convenient for
   (e.g.) highlighting a particular subplot by changing its color.

Note some of these settings can also be controlled using, e.g.,
``ax.format(title_kw={'weight':'bold'})`` instead of
``ax.format(rc_kw={'titleweight':'bold'})``.

To access a single setting, use ``rc.name`` or ``rc[name]``. To
access a group of setting by category name (e.g., the ``rcParams`` that
look like ``'axes.something'``), use ``rc.axes`` and a **dictionary**
will be returned.

To reset everything to the default state, use `~proplot.rcmod.rc_configurator.reset`. By
default, **settings are reset every time a figure is drawn** -- that is, when
a figure is rendered by the matplotlib backend or saved to file.

