General
=======

The subplots command
--------------------
Figures and axes with special ProPlot features can be created with `~proplot.subplots.subplots`. This function was inspired by `matplotlib.pyplot.subplots` and is packed with new features. It returns a figure instance and a special container "`~proplot.subplots.axes_list`" for the axes.

To generate complex grids, pass a 2D array of numbers corresponding to unique subplots. Use zero to allot empty space. For example, the following command creates a grid with one tall plot on the left, and two smaller plots on the right:

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots([[1,2],[1,3]])

while this command creates a grid with one long plot on top and two smaller plots on the bottom, with a gap in the middle (the zero):

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots([[1,1,1],[2,0,4]])

The corresponding axes "numbers" are added as attributes, and can be used for a-b-c labeling with the `~proplot.axes.BaseAxes.format` method.

For a single-axes figure (1 row, 1 column), simply call `~proplot.subplots.subplots` with no arguments. For simple grids, use the `nrows` and `ncols` keyword args; for example

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots(nrows=2)

The `~proplot.subplots.axes_list` container lets you invoke methods in every axes in ``axs`` simultaneously; for example:

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots()
   axs.format(xlim=(0,1)) # set to be identical everywhere
   axs[:3].format(color='red')  # make axis labels, spines, ticks red
   axs[3:].format(color='blue') # same, but make them blue

There are further arguments for drawing axes or figure "panels", setting particular axes
as map projections, sharing axis limits and tick/axis labels, scaling the
relative widths of rows/columns, setting axes aspect ratios and variable inter-axes
spacing, and setting average axes width/height instead of the full figure dimensions.

See `~proplot.subplots.subplots` for a full description of the arguments.

The subplot layout
------------------

If you specify just one of ``width``, ``height``, ``axwidth`` (fixes
figure width), or ``axheight`` (fixes figure height), the unspecified
dimension will be **scaled** such that the top-left subplot has aspect
ratio ``aspect``.

This is accomplished with the new `~proplot.gridspec.FlexibleGridSpec` class, subclassed from matplotlib’s `~matplotlib.gridspec.GridSpec` class. The “actual” ``wspace`` and ``hspace`` passed to `~matplotlib.gridspec.GridSpec` are zero – the spaces you see in your figure are empty subplot slots masquerading as spaces, whose widths are controlled by ``width_ratios`` and ``height_ratios``. Check out `FlexibleGridSpec.__getitem__`.

ProPlot also implements a new "tight layout" feature, heavily powered by `~proplot.gridspec.FlexibleGridSpec`. The `~proplot.subplots.Figure.smart_tight_layout` method is called whenever the figure is drawn (i.e. when it is rendered by the backend or saved to file).

The native matplotlib version, `~matplotlib.figure.Figure.tight_layout`, fits the figure borders over a box that perfectly encompasses all text, subplots, etc. However, there are three major problems with `~matplotlib.figure.Figure.tight_layout`:

1. It does not preserve subplot aspect ratios.
2. It cannot correct overlapping content or excessive whitespace inside the figure.
3. It can lead to overlapping content inside the subplot where there was none before, because the ``wspace`` and ``hspace`` coordinates are *relative* to the average subplot size while things like lines and text are specified in *absolute* units: "points" (i.e. 1/72 inches).

The new `~proplot.subplots.Figure.smart_tight_layout` method draws a tight bounding box that *preserves* panel widths and subplot aspect ratios. Not only that, it also adjusts the *spacing* between subplots, so there is no excessive whitespace and no overlap between axes tick labels and whatnot.

This is best demonstrated by example; see the :ref:`Tutorial` for more.

Axis sharing and spanning labels
--------------------------------
I've expanded the matplotlib `shared axis <https://matplotlib.org/examples/pylab_examples/shared_axis_demo.html>`_ setting and added a new "spanning axis label" feature. These are best demonstrated by example; see the :ref:`Tutorial` and `~proplot.subplots.subplots`.


Axes classes and the format method
----------------------------------
The `~proplot.subplots.subplots` method populates the `~proplot.subplots.Figure` object with either of three types of axes:

* `~proplot.axes.XYAxes`
* `~proplot.axes.CartopyAxes`
* `~proplot.axes.BasemapAxes`

Each of these inherits from the base class `~proplot.axes.BaseAxes`.

The most important new method you need to know is `~proplot.axes.BaseAxes.format`. This is your one-stop-shop for changing axis labels, tick labels, titles, etc. The keyword args passed to this function are interpreted as follows:

1. Any keyword arg matching the name of a custom ProPlot or native matplotlib "rc" settings will be applied to the axes (see the `~proplot.rcmod` documentation). If the name has "dots", simply omit them -- for example, ``title.weight`` becomes ``titleweight``.
2. Remaining keyword args are passed to the ``smart_update`` methods of the top-level class -- that is, the `~proplot.axes.XYAxes`, `~proplot.axes.CartopyAxes`, `~proplot.axes.BasemapAxes`. Use these to change settings specific to Cartesian axes or specific to map projections, like tick locations or toggling geographic features.
3. Finally, the remaining keyword args are passed to the `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.smart_update` method. This one controls "universal" settings -- namely, various titles and a-b-c label options.

Refer to the documentation on the different ``smart_update`` methods for usage information.

Some might argue that this method just replicates features already available from matplotlib -- so, some motivation is in order. To modify an axes property (e.g. an *x*-axis label) with the default API, you normally have to use a bunch of one-liner `~matplotlib.pyplot` commands, or method calls on axes and axis instances. This can get quite repetitive and quite verbose, resulting in lots of ugly, cumbersome boilerplate code.

Now, you just pass these settings to `~proplot.axes.BaseAxes.format`. Instead of having to remember the name of the function, whether it's attached to the `~matplotlib.pyplot` module or an object instance, and the order and names of the arguments, typing out a verbose command every time you want to change one little thing, you just pass a single keyword arg to `~proplot.axes.BaseAxes.format`.

Example:

.. code-block:: python

   import proplot as plot
   f, ax = plot.subplots()
   ax.format(xlabel='time (seconds)', ylabel='temperature (K)', title='20th century sea-surface temperature')


Note there is also the special `~proplot.axes.PanelAxes` class used for panels -- this class inherits from `~proplot.axes.XYAxes`, and only differs
in that (by default) calling `~proplot.axes.PanelAxes.legend` and
`~proplot.axes.PanelAxes.colorbar` on these axes will "fill" the entire axes with a legend or a colorbar (refer to the :ref:`Tutorial` to see this in action).


The rc object
-------------
A special object named `~proplot.rcmod.rc`, belonging to the
`~proplot.rcmod.rc_configurator` class, is created whenever you import ProPlot. This object gives you advanced control over the look of your plots.
**Use** `~proplot.rcmod.rc` **as your one-stop shop for changing global settings**.

The `~proplot.rcmod.rc` object controls built-in
`~matplotlib.rcParams` settings, a few custom :ref:`rcParams_new` settings,
and some magic :ref:`rcGlobals` settings that apply to groups of other
settings and keep them synced. Tables of these settings are found in the `~proplot.rcmod` documentation.

To modify any :ref:`rcGlobals`, :ref:`rcParams_new`, or `~matplotlib.rcParams` setting, you have four options:

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

"Outer" Panels
--------------

ProPlot figures may have “outer” panels on the bottom, left,
or right of the figure, accessed with ``fig.bottompanel``,
``fig.leftpanel``, and ``fig.rightpanel``, respectively. They are instances of `~proplot.axes.PanelAxes`, and can be set up with `~proplot.subplots.subplots` in several ways:

-  ``bottompanel=True``: Allot space for a single panel spanning all
   columns of subplots.
-  ``bottompanels=True``: Allot space for ``n`` separate panels occupying
   the ``n`` columns of subplots.
-  ``bottompanels=[array]``: Allot space for an arbitrary number of
   panels spanning contiguous subplot columns – for example,
   ``bottompanels=[1,1,1,2,2,3]`` draws 3 panels, spanning 3 columns, 2
   columns, and 1 column respectively.

To access the nth panel, use ``fig.bottompanel[n]``. If you didn't request a panel in the call to `~proplot.subplots.subplots`, this will raise a ``"Panel does not exist."`` attribute error.

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

Outer panel settings can be configured with a bunch of `~proplot.subplots.subplots` keyword arguments; see the :ref:`Documentation` for details.

"Inner" Panels
--------------

ProPlot also provides utilities for making “**inner panels**”. These may
be useful where you want a colorbar for every plot, a legend outside of
every axes, or want to show the x/y-direction statistics for some 2D
value plotted in your subplot (e.g. the x-direction mean, variance,
etc.).

The procedure for requesting inner panels is similar:

* ``innerpanels='r'`` draws panels on the right of each subplot, ``innerpanels='rt'`` draws panels on the top and the right.
* ``innerpanels={1:'r', (2,3):''}`` or ``innerpanels={range(5):'bt', 5:''}`` draws inner panels for *particular subplot numbers*.

As with the outer panels, you can also use ``innercolorbars='r'`` to
draw panels with default widths/spacing suitable for colorbars.

Inner panel settings can be configured with a bunch of `~proplot.subplots.subplots` keyword arguments; see the :ref:`documentation` for details.

Colorbar Enhancements
---------------------

See `~proplot.axes.colorbar_factory` for more info. Normally, to create a colorbar, you need a “mappable” instance – i.e.,
something with a `get_cmap` method, returned by `contourf`,
`pcolor`, etc.

With ProPlot, you have the following two additional options:

* ``colorbar(colors, values=values)`` creates a `~matplotlib.colors.ListedColormap` using a list of color strings or ``(R,G,B)`` tuples, then draws a colorbar with the values in the iterable `values` mapped to each color.
* ``colorbar(handles, values=values)`` infers colors from a list of “plot handles” -- i.e. anything with a `get_color` method, including handles returned by `~matplotlib.axes.Axes.plot`).

Two additional options make it easy to configure your colorbar geometry:

* ``length=fraction``, where ``0 <= fraction <= 1``, will make the colorbar span a *fraction* of the horizontal/vertical extent of the axes it is filling.
* ``extendlength=size`` controls the length of the “triangles” representing out-of-bounds colors (drawn when you use ``extend='min'``, ``extend='max'``, or ``extend='both'``). Since the “triangles” are now specified in physical units, they will always match other colorbars in the figure.

Legend Enhancements
-------------------

See `~proplot.axes.legend_factory` for more info. Normally, the legend handles are plotted in column-major order by default. Now you can choose between row-major and column-major, and the
former is the new default.

You can also create pseudo-legends with handles that **are not** aligned by column, using
``align=False`` -- or by passing a *list of lists of handles* instead of a list of handles. This actually creates a bunch of centered single-row legends stacked
on top of each other. This can be handy when you want to organize rows of handles logically,
or you have just a couple handles beneath a really long row and you want them centered.


Cartesian and map axes
======================

Formatting Cartesian axes
-------------------------
The following is a brief overview of valid arguments for the `~proplot.axes.XYAxes`
`~proplot.axes.XYAxes.smart_update` method (which receives arguments from the `~proplot.axes.BaseAxes.format` command). The below just highlights some of the more useful ones; refer to the documentation for a complete description.

You can very quickly modify tick positions with the `xlocator` and `ylocator` keywords (or their aliases, `xticks` and `yticks`). These accept a number of possible arguments:

*  A number (e.g. ``xticks=N``) ticks every N data values.
*  A string will look up any of the `matplotlib.ticker`
   locators by key name, e.g. ``'log'``.
*  A list of numbers will tick those specific locations.

I recommend using `plot.arange` to generate lists of ticks –
it’s like `numpy.arange`, but is **endpoint-inclusive**, which more often than
not is what you'll want in this context.

You can also control the tick label format with `xformatter` and `yformatter` keywords (or their aliases, `xticklabels` and `yticklabels`). These accept a number of possible arguments:

* ``'%.f'`` for classic `%-style formatting <https://pyformat.info/>`_, or ``{}`` for newer `'string'.format(value)` formatting.
* ``'lat'``, ``'deglat'``, ``'lon'``, ``'deglon'``, and ``'deg'``
  format axis labels with cardinal direction indicators and/or degree
  symbols, as denoted by the respective names.
* ``'pi'``, ``'e'``, ``('symbol', scale)`` will format tick labels represented as
  fractions of some symbol (the first 2 are :math:`\pi` and Euler's constant, provided for convenience).
* A list of strings (e.g. ``xticklabels=['a', 'b', 'c']``) will simply label existing ticks with that list.

You can control the axis scale (e.g. ``'linear'`` vs. ``'log'``) with the `xscale` and `yscale` keywords. There are also some new scales available, described below:

-  The "inverse" scale ``'inverse'``. Useful for, e.g., having
   wavenumber and wavelength on opposite sides of the same plot.
-  The sine-weighted and "Mercator" axis scales, ``'sine'`` and
   ``'mercator'``.
-  The "cutoff" scale, allowing arbitrary
   zoom-ins and zoom-outs over segments of an axis. This is actually an
   axis scale created **on-the-fly**, invoked with the tuple
   ``('cutoff', scale, lower, upper)`` where ``lower``
   and ``upper`` are the boundaries within which the axis scaling is
   multiplied by ``scale``. Use ``np.inf`` for a hard cutoff, or
   use ``('cutoff', scale, position)`` to scale every coordinate after
   position ``position`` by ``scale``.


Map projections
---------------
ProPlot also lets you set up axes with geographic projections using either of 2 packages: `~mpl_toolkits.basemap.Basemap` or `~cartopy.crs.Projection`.
See the `~proplot.projs` documentation for a handy table of available projection names.

Note that `~mpl_toolkits.basemap` is `no longer under active development <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`_ -- cartopy is the intended replacement, as it is integrated more intelligently with the matplotlib API.
However, for the time being, basemap retains one advantage over cartopy. Namely, `support for labeling meridians and parallels <https://github.com/SciTools/cartopy/issues/881>`_. I therefore decided to support both, for the time being.

Projections are configured with the ``proj`` and ``proj_kw`` keyword args via the `~proplot.subplots.subplots` command. Set the map projection for all subplots with ``proj='proj'``, or separately for different subplots with e.g. ``proj={1:'proj1', (2,3):'proj2', 4:'name3'}``. In the latter case, the integers and integer tuples correspond to **the subplot number**.

In the same way, you can pass keyword args to the cartopy `~cartopy.crs.Projection` and `~mpl_toolkits.basemap.Basemap` class initializers using ``proj_kw={'name':value}`` or e.g. ``proj_kw={1:'proj1', (2,3):'proj2'}``.

You can also choose between cartopy and basemap using ``basemap=False`` or e.g. ``basemap={1:True, 2:False}``.

As a simple example, the following creates 3 side-by-side `Hammer projections <https://en.wikipedia.org/wiki/Hammer_projection>`_ using Cartopy.

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots(ncols=3, proj='hammer', basemap=False)

Cartopy axes
------------
When you specify the ``proj`` keyword arg with ``basemap=False``, a `~proplot.axes.CartopyAxes` instance (subclassed from the cartopy `~cartopy.mpl.geoaxes.GeoAxes` class) is created. As shown above, you can now declare the projection by **string name**, instead of having to reference cartopy `~cartopy.crs.Projection` classes directly.

In cartopy, you usually need to supply ``transform=cartopy.crs.PlateCarree()`` to the plotting method (see `this example <https://scitools.org.uk/cartopy/docs/v0.5/matplotlib/introductory_examples/03.contours.html>`_). With ProPlot, **this is done by default**.

Other aspects of cartopy axes can be controlled with `~proplot.axes.CartopyAxes.smart_update` (i.e. with `~proplot.axes.BaseAxes.format`, which calls `~proplot.axes.CartopyAxes.smart_update`).

Basemap axes
------------
When you specify the ``proj`` keyword arg with ``basemap=True``, a `~proplot.axes.BasemapAxes` instance is created. This class allows you to access basemap plotting utilities **directly on the axes as a method**, instead of having to call the method from the `~mpl_toolkits.basemap.Basemap` instance.

To fix issues with the "seam" on the edge of the map, I've overridden several plotting methods on `~proplot.axes.BasemapAxes` -- data will be automatically circularly rolled until the left-hand-side comes after the map seam, then interpolated to the seam longitudes.

Other aspects of basemap axes can be controlled with `~proplot.axes.BasemapAxes.smart_update` (i.e. with `~proplot.axes.BaseAxes.format`, which calls `~proplot.axes.BasemapAxes.smart_update`).

Colormaps and colors
====================

A figure prepared for publication should be a work of art. Your
figures should tell the entire story – the article text just fills in the blanks.
Several tools have been added to help make your graphics both visually
appealing and informative.

Colormaps
---------

A **new colormap class** analogous to `~matplotlib.colors.LinearSegmentedColormap` is now
available, called `~proplot.colortools.PerceptuallyUniformColormap`. This class linearly
interpolates through hue, chroma, and luminance space (with hues allowed
to vary circularly), instead of RGB space as with
`~matplotlib.colors.LinearSegmentedColormap` 

The colors in a `~proplot.colortools.PerceptuallyUniformColormap` can span either of `4
HSV-like colorspaces <http://www.hsluv.org/comparison/>`__: classic HSV,
perceptually uniform HCL, or HSLuv/HPLuv (which are forms of HCL adapted
for this kind of usage).

I’ve packaged several `~proplot.colortools.PerceptuallyUniformColormap` maps
with ProPlot, along with perceptually uniform maps from several other projects.
:ref:`Colormaps` provides a table of the registered colormaps.

Colormap generation
-------------------

Generate a `~proplot.colortools.PerceptuallyUniformColormap` on-the-fly by passing a
**dictionary** to any plotting function that accepts the `cmap`
keyword argument.

Example:

.. code:: python

   f, ax = plot.subplots()
   ax.contourf(np.random.rand(10,10), cmap={'h':['red', 'red+30'], 'c':80, 'l':[50, 100], 'space':'hpl'}

The arguments can be single numbers, lists of numbers, or single/lists
of **color strings**. In the latter case, the corresponding channel
value (hue, chroma, or luminance) for that color will be looked up and
applied. You can end any color string with ``+N`` or ``-N`` to offset
the channel value by the number ``N``, as shown above. Note you can also
use ``s`` (saturation) instead of ``c`` (chroma), or the full names
``'hue'``, ``'chroma'``, ``'luminance'``, or ``'saturation'``.

Create single-hue colormaps on-the-fly by passing a string that looks
like ``cmap='name'`` or ``cmap='nameXX'``, where ``name`` is any
registered color string (the corresponding hue will be looked up) and
``XX`` is the lightness value.

Example:

.. code:: python

   f, ax = plot.subplots()
   ax.pcolormesh(np.random.rand(10,10), cmap='sky blue70', cmap_kw={'name':'my_cmap', 'save':True})

creates a monochrome colormap. It also saves the colormap with the name
``'my_cmap'``, using the `cmap_kw` dictionary argument.

The default colormap can be set with ``plot.rc.cmap = <cmap spec>`` or ``plot.rc.cmap = (<cmap spec>, <cmap kwargs>)``,
where the colormap specification (and optional keyword args) are passed through `~proplot.Colormap`.

Color cycles
------------

In addition to the new colormaps, new “color cycles” are also available
(i.e. the automatic color order used for drawing multiple lines).
:ref:`Color cycles` provides a table of these cycles.

The default cycler can be set with ``plot.rc.cycle = <cycle spec>`` or ``plot.rc.cycle = (<cycle spec>, <cycle kwargs>)``,
where the cycle specification (and optional keyword args) are passed through `~proplot.Cycle`.
The cycler can also be temporarily changed by passing ``cycle='name'`` (and, optionally, ``cycle_kw={'key':value}``)
to any plotting command that ordinarily loops through a color cycle, e.g. ``plot`` and ``bar``.


The **distinction between a “colormap” and “color cycle” is now fluid**:

1. All color cycles are implemented as `~matplotlib.colors.ListedColormap` instances; you can request them as colormaps with ``cmap='cycle_name'``.
2. Cycles can be generated on the fly from the colormaps by specifying e.g. ``cycle=('cmap_name', N)``, where ``N`` is the number of colors over the registered colormap you’d like to sample. If you just use ``cycle='cmap_name'``, the default is 10 colors.


The following generates a cycle of 5 colors over the matplotlib builtin colormap
``'blues'``, excluding the very brightest colors:

.. code:: python

   f, ax = plot.subplots()
   ax.plot(np.random.rand(10,5), cycle=('blues', 5), cycle_kw={'x':(0.2,1)})


Registered color names
----------------------

New colors names have been added from the `XKCD color-naming
project <https://xkcd.com/color/rgb/>`__, so-called “crayon” colors
provided with `Seaborn <https://seaborn.pydata.org/>`__, and Open Color
web-design color palette. Colors that aren't sufficiently perceptually
distinct are eliminated, so it's easier to pick from the color table.

:ref:`Table of colors` provides table of the newly registered colors.

Contour and pcolor
------------------

This one is a small change – I've fixed the well-documented `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__ and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__ issues by automatically changing the edge colors when `contourf`, `pcolor`, and `pcolormesh` are called.
