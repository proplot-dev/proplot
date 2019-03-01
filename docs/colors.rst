Color management
================

A figure prepared for publication should be a work of art. Your
figures should tell the entire story – the article text just fills in the blanks.
Several tools have been added to help make your graphics both visually
appealing and informative.

Colormaps
---------

A **new colormap class** analagous to `~matplotlib.colors.LinearSegmentedColormap` is now
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
web-design color palette. Colors that aren’t sufficiently perceptually
distinct are eliminated, so it's easier to pick from the color table.

:ref:`Color names` provides table of the newly registered colors.

Contour and pcolor
------------------

This one is a small change – I’ve fixed the well-documented
`white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__
and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__
issues by automatically changing the edgecolors when `contourf`,
`pcolor`, and `pcolormesh` are called.
