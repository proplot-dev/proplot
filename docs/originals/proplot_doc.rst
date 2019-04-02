{% comment %} \* TOC {% endcomment %} {% comment %} {:toc} {% endcomment
%}

For worked examples of the features described below, see the
`showcase <%7B%7B%20site.baseurl%20%7D%7D%7B%%20link%20_tools/proplot.md%20%%7D>`__.
Please note that some features have not been covered yet; check back
later.

Import as follows:

{% comment %} I recommend importing with {% endcomment %}

.. code:: python

   import proplot as plot

Figure creation
===============

The ``plot.subplots`` command was inspired by the ``matplotlib.pyplot``
command of the same name. This version is packed with new features.

The subplots command
--------------------

To generate complex grids, pass a 2D array of numbers corresponding to
unique subplots. Use zero to allot empty space.

Example:

.. code:: python

   f, axs = plot.subplots([[1,2],[1,3]])

creates a grid with one tall plot on the left, and two smaller plots on
the right, while

.. code:: python

   f, axs = plot.subplots([[1,1,1],[2,0,4]])

creates a grid with one long plot on top and two smaller plots on the
bottom with a gap in the middle (the zero). The corresponding
``number``\ s will be added as axes attributes; you can use ``format``
to trigger a-b-c labeling based on these numbers (see below).

Pass no arguments to generate default single-axes figure (1 row, 1
column), or use ``nrows`` and ``ncols`` to generate simple subplot
grids.

Example:

.. code:: python

   f, axs = plot.subplots(nrows=2)

creates two rows with one column. Note that if you did not use
``array``, the subplots are automatically numbered in row major order.

The two return values are the figure instance, ``f``, and a list of
axes, ``axs``. This list of axes **is no ordinary list** – it is a
special class called ``axes_list`` that allows you to bulk-invoke
methods on multiple axes simultaneously, no for-loop necessary.

Example:

.. code:: python

   f, axs = plot.subplots()
   axs.format(xlim=(0,1)) # set to be identical everywhere
   axs[:3].format(color='red') # make axis labels, spines, ticks red
   axs[3:].format(color='blue') # same, but make them blue

Arguments
---------

The ``subplots`` command takes a bunch of handy arguments, described in
detail below. Note that, while *numeric* sizing arguments are assumed to
be in inches, you can specify sizes with **arbitrary units** using e.g.
``width='12cm'``, ``wspace='5mm'``, ``hspace='3em'`` (3 em squares),
etc.

-  ``width``, ``height``: Sets the figure width and height.
-  ``axwidth``, ``axheight``: Sets the average width, height of your
   axes. This is convenient where you don’t care about the figure
   dimensions, but just want your axes to have enough “room”. Default is
   ``axwidth=2``.
-  ``aspect``: The aspect ratio (width over height) of the **top-left**
   subplot region, specified using a number or length-2 iterable (e.g.
   ``aspect=1.5`` or ``aspect=(3,2)``). Default is ``1``
-  ``wspace``, ``hspace``: The width/height spacing between axes.
-  ``wratios``, ``hratios``: The width/height ratios for columns/rows of
   axes.
-  ``top``, ``bottom``, ``left``, ``right``: The border space. If you
   have outer panels (see below), these apply to the space between outer
   panels and the main subplot region.

Journal standards
-----------------

To create figures with dimensions that satisfy journal standards, use
the ``journal`` keyword argument.

Example:

.. code:: python

   f, axs = plot.subplots(ncols=3, nrows=2, journal='ams2') # medium-sized figure for AMS journal

The currently available specifiers are as follows; feel free to contact
me with additional standards, and I will add them. \* ``pnas1``,
``pnas2``, ``pnas3``: Single-column, medium, and two-column figure
widths, from `here <http://www.pnas.org/page/authors/submission>`__. \*
``ams1``, ``ams2``, ``ams3``, ``ams4``: Single-column, medium,
two-column, and landscape-page figure widths, from
`here <https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/>`__.
\* ``agu1``, ``agu2``, ``agu3``, ``agu4``: Single-column, medium,
two-column, and landscape-page figure sizes, from
`here <https://publications.agu.org/author-resource-center/figures-faq/>`__.

Smarter subplot layout
----------------------

If you specify just one of ``width``, ``height``, ``axwidth`` (fixes
figure width), or ``axheight`` (fixes figure height), the unspecified
dimension will be **scaled** such that the top-left subplot has aspect
ratio ``aspect``. **Spacing and panel widths are held fixed** during
this scaling.

This is accomplished with the new ``FlexibleGridSpec`` class, subclassed
from matplotlib’s ``GridSpec`` class. The “actual” ``wspace`` and
``hspace`` passed to ``GridSpec`` are zero – the spaces you see in your
figure are empty subplot slots *masquerading* as spaces, whose widths
are controlled by ``width_ratios`` and ``height_ratios``. Check out
``FlexibleGridSpec.__getitem__``.

Note this also means **inter-subplot spacing is now variable**. You can
specify ``wspace`` and ``hspace`` with: 1. A scalar constant, e.g.
``wspace=0.2``. 1. A list of different spacings, e.g.
``wspace=[0.1,0.5]`` to offset the 3rd column in a 3-column plot from
the rest.

Smarter “tight” layout
----------------------

The subplot layout changes allowed me to create the
``smart_tight_layout`` method. By default, this method is **called
whenever the figure is drawn** (i.e. when it is rendered by the
matplotlib backend or saved to file) – disable this behavior with
``plot.subplots(tight=False)``.

Previously, ``tight_layout`` could be used to fit the figure borders
over a box that perfectly encompasses all artists (i.e. text, subplots,
etc.). However, because ``GridSpec`` spaces are relative to the subplot
dimensions, changing the figure dimensions also changes the
inter-subplot spacings. Since your font size is specified in points
(i.e. a physical unit), *this can easily cause text to overlap with
other subplots where it didn’t before*.

The new ``smart_tight_layout`` method draws a tight bounding box that
**preserves inter-subplot spacing, panel widths, and subplot aspect
ratios**.

Panels
======

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

Cartopy + Basemap integration
=============================

For projection subplots, specify ``projection='name'`` with either
``package='basemap'`` or ``package='cartopy'``.

Extra arguments to ``subplot`` will be passed to the ``basemap.Basemap``
and ``cartopy.crs.Projection`` classes (the relevant cartopy class will
be selected based on the ``'name'`` string).

Control the map projection type with ``proj='proj'`` or e.g.
``proj={1:'proj1', (2,3):'proj2', 4:'name3'}``. In the latter case, the
integers and integer tuples correspond to **the subplot number**.

In the same way, you can pass keyword arguments (e.g. ``lon_0``) to the
``cartopy.crs.Projection`` or ``basemap.Basemap`` classes using
``proj_kw={'name':value}`` or e.g.
``proj_kw={1:'proj1', (2,3):'proj2'}``. You can also choose between
cartopy and basemap using ``basemap=False`` or e.g.
``basemap={1:True, 2:False}``.

Example:

::

   f, ax = plot.subplots(ncols=3, proj='hammer', basemap=False)

creates 3 side-by-side `Hammer
projections <https://en.wikipedia.org/wiki/Hammer_projection>`__ using
Cartopy.

Note that `Basemap is no longer under active
development <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__
– cartopy is integrated more intelligently with the matplotlib API.
However, for the time being, basemap *retains several advantages* over
cartopy (namely `more tools for labeling
meridians/parallels <https://github.com/SciTools/cartopy/issues/881>`__
and more available projections – see
`basemap <https://matplotlib.org/basemap/users/mapsetup.html>`__ vs.
`cartopy <https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html>`__).
Therefore, I decided to support both.

Cartopy axes
------------

When you use ``proj=<something>`` with ``basemap=False``, a
``plot.CartopyAxes`` instance (subclassed from the cartopy ``Geoaxes``
class) is created. As shown above, you can now declare the projection by
**string name**, instead of having to reference
``cartopy.crs.ProjectionName`` classes directly.

In cartopy, you usually need to `supply
``transform=cartopy.crs.PlateCarree()`` to the plotting
method <https://scitools.org.uk/cartopy/docs/v0.5/matplotlib/introductory_examples/03.contours.html>`__.
With ProPlot, this is done by default.

Other aspects of cartopy axes can be controlled with ``ax.format()`` –
see the “Formatting” section.

Basemap axes
------------

When you use ``proj=<something>`` with ``basemap=True``, a
``plot.BasemapAxes`` instance is created. This class allows you to
access basemap plotting utilities **directly on the axes as a method**,
instead of having to call the method from the ``Basemap()`` instance.

To fix issues with the “seam” on the edge of the map, I’ve overridden
several plotting methods on ``BasemapAxes`` – data will be automatically
circularly rolled until the left-hand-side comes after the map seam,
then interpolated to the seam longitudes.

Other aspects of basemap axes can be controlled with ``ax.format()`` –
see the “Formatting” section.

Formatting
==========

The new ``format`` method, available on every axes returned from
``subplots()``, is a versatile and powerful tool. This section describes
its uses.

Motivation
----------

To modify an axes property (e.g. an x-axis label) with the default API,
you normally have to use a bunch of one-liner ``pyplot`` commands (or
method calls on axes/axis objects). This can get repetitive and quite
verbose, resulting in lots of boilerplate code.

Now, you can pass these settings to ``format``. Instead of having to
remember the name of the function, whether it’s attached to ``pyplot``
or an object instance, and the order/names of the arguments, you just
have to remember one thing – the name of the keyword argument.

The ``format`` method also abstracts away some inconsistencies and
redundancies in the matplotlib API – now, *There’s Only One (obvious)
Way To Do It*.

Example:

.. code:: python

   ax.format(xlabel='time (seconds)', ylabel='temperature (K)', title='20th century sea-surface temperature')

.. _arguments-1:

Arguments
---------

The ``format`` method has a ton of optional arguments, described in
detail below. Note wherever you see ``_kw``, this indicates an optional
dictionary that will be supplied as keyword arguments to some function –
``ax.text`` for labels and stuff, ``plot.locator`` for tick locators,
``plot.formatter`` for tick formatters, and ``plot.scale`` for axis
scales.

Note the same convention is used in other situations: e.g., ``cmap_kw``
is passed to ``plot.colormap`` when passed to a command that accepts
``cmap``, and ``cycle_kw`` is passed to ``plot.cycle`` (alias
``plot.colors``) when passed to a command that accepts ``cycle``.

**Titling arguments**: \* ``suptitle``, ``suptitle_kw``: The “super”
title (axes-spanning title). \* ``title``, ``titlepos``, ``title_kw``:
The axes title and title position, specified with an optional string up
to 2 characters wide, where ``i`` denotes a title *inside* axes bounds,
``o`` *outside* (the default), and ``l``, ``c``, and ``r`` denote titles
aligned to the left, center (default), and right side of the axes. \*
``abc``, ``abcpos``, ``abcformat``, ``abc_kw``: The “a-b-c” subplot
labeling and label position (specified in the same way as the title
position). Use e.g. ``abcformat='(a)'``, ``abcformat='a.'``, etc. to
format the labels

**Axis-related arguments**: \* ``xlabel``, ``ylabel``, ``xlabel_kw``,
``ylabel_kw``: These control the axis labels and formatting. \*
``xlim``, ``ylim``, ``xreverse``, ``yreverse``: Use these to control the
axis data limits. Use ``xreverse`` and ``yreverse`` to reverse the axis
limits. \* ``xtickrange``, ``ytickrange``: Use these to restrict the
range over which your major ticks are labeled (e.g.,
``xtickrange=(-1,1), xlim=(-5,5)`` if you only want enough labels to
denote the major tick step size). \* ``xscale``, ``yscale``,
``xscale_kw``, ``yscale_kw``: These control the axis scales – use e.g.
``log``, ``linear``, ``('cutoff', 10, 5)`` (see below for details). The
keyword arguments get passed to the scale constructor. \* ``xspineloc``,
``yspineloc``, ``xloc``, ``yloc``: One of ``bottom`` or ``top`` (for
``xspineloc``), ``left`` or ``right`` (for ``yspineloc``), ``both``, or
``neither`` – controls where we draw axes spines. The ``xloc`` and
``yloc`` arguments are just aliases \* ``xgrid``, ``ygrid``, ``grid``:
True or False – whether to draw grid lines. Use ``grid`` to set both the
``x`` and ``y`` properties. \* ``xtickloc``, ``ytickloc``: One of
``bottom`` or ``top`` (for ``xtickloc``), ``left`` or ``right`` (for
``ytickloc``), ``both``, or ``neither``. \* ``xtickdir``, ``ytickdir``,
``tickdir``: One of ``in``, ``out``, or ``inout`` – controls which
direction ticks point on the spine. \* ``xtickminor``, ``ytickminor``,
``tickminor``, ``xgridminor``, ``ygridminor``, ``gridminor``: True or
False – whether to draw minor ticks (default True) and minor tick grid
lines (default False). \* ``xticklabeldir``, ``yticklabeldir``: One of
``in`` or ``out`` – whether to place tick label text inside or outside
the axes. \* ``xlocator``, ``xminorlocator``, ``ylocator``,
``yminorlocator``, ``xticks``, ``yticks``, ``xminorticks``,
``yminorticks``: These are flexible arguments for setting up major and
minor tick positions (the ``ticks`` arguments are aliases for the
``locator`` arguments). See below for details. Use ``xlocator_kw``,
``xminorlocator_kw``, ``ylocator_kw``, ``yminorlocator_kw`` to pass
extra arguments to the locator constructor. \* ``xformatter``,
``yformatter``, ``xticklabels``, ``yticklabels``: These are flexible
arguments for setting up major tick labels (the ``ticklabels`` arguments
are aliases for the ``formatter`` arguments). See below for details. Use
``xformatter_kw``, ``yformatter_kw`` to pass extra arguments to the
formatter constructor. \* ``latlabels``, ``lonlabels``, ``ylabels``,
``xlabels`` (maps only): These control which sides of your map
projection that meridian and parallel labels will appear. Use e.g.
``latlabels='lb'`` to label parallels on the bottom and left axes. The
latter two are simple aliases. \* ``latlocator``, ``latminorlocator``,
``lonlocator``, ``lonminorlocator`` (maps only): Control major and minor
grid lines. You may also use ‘x’ instead of ‘lon’ and ‘y’ instead of
‘lat’. \* ``land``, ``ocean``, ``land_kw``, ``ocean_kw`` (maps only):
Whether to fill in land masses or oceans with a patch. The patch
properties can be controlled with the ``kw`` dictionaries or using the
``rc`` configurator (see below).

The ``format()`` method also handles ``rc`` configuration; see the
“Settings Management” section for more information.

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
I’ve also created several special formatter classes: \* Basic formatters
``'%.f'``, for classic ```%``-style
formatting <https://pyformat.info/>`__, or ``{x}`` for newer
``'string'.format(x=value)`` formatting. \* Coordinate formatters
``'lat'``, ``'deglat'``, ``'lon'``, ``'deglon'``, ``'deg'``: Use these
to format axis labels with cardinal direction indicators and/or degree
symbols, as denoted by the respective names. \* Fraction formatters
``'pi'``, ``'e'``, ``('symbol', scale)``: For tick labels represented as
fractions of some symbol. \* Explicit tick labels with e.g.
``xticklabels=['a', 'b', 'c']`` – adds specific strings to existing
major ticks.

Advanced settings management
============================

A special object named ``rc`` (belonging to a class called
``rc_configurator``) is created whenever you import ProPlot. This object
gives you advanced control over the look of your plots. **Use**
``plot.rc`` **as your one-stop shop for changing global settings**. ##
The ``rc`` object The ``rc`` object can be used to change built-in
``matplotlib.rcParams`` settings, a few custom “``rcSpecial``” settings,
and some magic “``rcGlobal``” settings that apply to groups of other
settings and keep them synced – e.g., tick, spine, and tick label
colors.

To modify ``rc`` settings, you have three options: 1. Change one global
setting using ``plot.rc.name = value`` or ``plot.rc['name'] = value``.
Note that, for ``rcParams`` settings with ‘dots’ in their name, you will
have to use ``plot.rc['category.name'] = value``. 1. Update several
global settings at once using
``plot.rc.update({'name1':value1, 'name2':value2})`` or
``plot.rc.update(name1=value1, name2=value2)``, just like you would
update a dictionary. 1. Change local settings using
``ax.format(rc_kw={'name1':value1, 'name2':value2})`` or
``ax.format(name1=value1, name2=value2)``. Note that, for this last
option, **the rc settings will only be applied locally** (i.e. to the
axes on which ``format()`` is being invoked). This can be convenient for
(e.g.) highlighting a particular subplot by changing its color.

To **access** a single setting, use ``rc.name`` or ``rc[name]``. To
access a group of setting by category name (e.g., the ``rcParams`` that
look like ``'axes.something'``), use ``rc.axes`` and a **dictionary**
will be returned.

Global settings
---------------

The following is an overview of the available “``rcGlobal``” settings.
\* ``color``, ``xcolor``, ``ycolor``: The color of axis spines, tick
marks, tick labels, and labels. Use the ``x`` and ``y`` versions to just
change settings for one axis. \* ``facecolor``, ``facehatch``: The
background color and optional hatching pattern (e.g. ``'xxx'``; see
`this
demo <https://matplotlib.org/gallery/shapes_and_collections/hatch_demo.html>`__)
for an axes. The latter is useful where you wish to highlight invalid
(transparent) “NaN” data in a ``pcolormesh`` or ``contourf`` plot. \*
``small``, ``large``: Font size for legend text, tick labels, axis
labels, and manually placed text created by ``ax.text`` (``small``), and
for titles, super titles, and a-b-c labels (``large``). \* ``fontname``:
The font name used for all text in the figure. ProPlot comes packaged
with a bunch of desirable fonts, and **changes the default from DejaVu
Sans (or Bitstream Vera) to Helvetica Neue**. When you first import
ProPlot, run ``plot.install_fonts()`` and restart your ipython session
to install them. \* ``linewidth``, ``minorwidth``: thickness of axes
spines and major tick lines (``linewidth``), and minor tick lines
(``minorwidth``, where minor tick line thickness =
``linewidth * minorwidth``). \* ``gridwidth``, ``gridratio``: thickness
of gridlines (``gridwidth``), and minor gridlines (``gridratio``, where
minor gridline thickness = ``gridwidth * gridratio``). \* ``gridalpha``,
``gridcolor``, ``gridstyle``: the transparency, color, and line style
for your major and minor gridlines. \* ``ticklen``, ``tickratio``:
length of major ticks (``ticklen``), and minor ticks (``tickratio``,
where minor tick lengths = ``ticklen * tickratio``). \* ``tickdir``:
tick direction, one of ``out``, ``in``, or ``inout``. \* ``abcweight``,
``titleweight``, ``suptitleweight``: the font weight (one of
``ultralight``, ``light``, ``normal``, ``medium``, ``demi``, ``bold``,
``very bold``, or ``black``) for your title text. Note that many fonts
only have ``normal`` or ``bold`` available; if you request another
weight, the “closest” availble weight will be selected.

Note some of these settings can also be controlled using, e.g.,
``ax.format(title_kw={'weight':'bold'})`` instead of
``ax.format(rc_kw={'titleweight':'bold'})``.

Use ``plot.rc.reset()`` to reset everything to the initial state. By
default, **settings are reset every time a figure is drawn** (i.e. when
a figure is rendered by the matplotlib backend or saved to file).

Colors and stuff
================

IMHO a figure prepared for publication should be a work of art. Your
figures tell the entire story – article text just fills in the blanks.
Several tools have been added to help make your graphics both visually
appealing and informative.

Colormaps
---------

A **new colormap class** analagous to ``LinearSegmentedColormap`` is now
available, called ``PerceptuallyUniformColormap``. This class linearly
interpolates through hue, chroma, and luminance space (with hues allowed
to vary circularly), instead of RGB space as with
``LinearSegmentedColormap``.

I’ve made several ``PerceptuallyUniformColormap``\ s already that come
packaged with ProPlot. I’ve also added the perceptually uniform
`cmOcean <https://matplotlib.org/cmocean/>`__ colormaps. **A table of
the new colormaps can be found at the end of
the**\ `showcase <%7B%7B%20site.baseurl%20%7D%7D%7B%%20link%20_tools/proplot.md%20%%7D>`__\ **.**

The colors in a ``PerceptuallyUniformColormap`` can span either of `4
HSV-like colorspaces <http://www.hsluv.org/comparison/>`__: classic HSV,
perceptually uniform HCL, or HSLuv/HPLuv (which are forms of HCL adapted
for this kind of usage). They can be specified using ``space='string'``
where ``'string'`` is any of ``'hsv'``, ``'hcl'``, ``'hsl'``, or
``'hpl'``.

Colormap generation
-------------------

Generate a ``PerceptuallyUniformColormap`` on-the-fly by passing a
**dictionary** to any plotting function that accepts the ``cmap``
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
``'my_cmap'``, using the ``cmap_kw`` dictionary argument.

Color cycles
------------

In addition to the new colormaps, new “color cycles” are also available
(i.e. the automatic color order used for drawing multiple lines). **A
table of the new color cycles can be found at the end of
the**\ `showcase <%7B%7B%20site.baseurl%20%7D%7D%7B%%20link%20_tools/proplot.md%20%%7D>`__\ **.**

The color cycler can be set with ``plot.rc.cycle = 'name'`` or by
passing ``cycle='name'`` to any command that plots lines/patches
(``plot``, ``bar``, etc.).

The **distinction between a “colormap” and “color cycle” is now fluid**:
0. All color cycles are defined as ``ListedColormap``\ s, and you can
request them as colormaps with ``cmap='cycle_name'``. 1. Cycles can be
generated on the fly from ``LinearSegmentedColormap``\ s by specifying
e.g. ``cycle=('cmap_name', N)`` where ``N`` is the number of colors over
the registered colormap you’d like to sample. If you just use
``cycle='cmap_name'``, the default will be 10 colors.

Example:

.. code:: python

   f, ax = plot.subplots()
   ax.plot(np.random.rand(10,5), cycle=('blues', 5), cycle_kw={'x':(0.2,1)})

generates a cycle of 5 colors over the matplotlib builtin colormap
``blues``, excluding the brightest/whitest colors (using the
``cycle_kw`` dictionary argument).

Registered color names
----------------------

New colors names have been added from the `XKCD color-naming
project <https://xkcd.com/color/rgb/>`__, so-called “crayon” colors
provided with `Seaborn <https://seaborn.pydata.org/>`__, and Open Color
web-design color palette. Colors that aren’t sufficiently perceptually
distinct are eliminated, so its easier to pick from the color table. **A
table of the new registered colors can be found at the end of
the**\ `showcase <%7B%7B%20site.baseurl%20%7D%7D%7B%%20link%20_tools/proplot.md%20%%7D>`__\ **.**

Contour and pcolor
------------------

This one is a small change – I’ve fixed the well-documented
`white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__
and
`white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__
issues by automatically changing the edgecolors when ``contourf``,
``pcolor``, and ``pcolormesh`` are called.
