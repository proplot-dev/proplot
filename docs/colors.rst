Colors and stuff
================

A figure prepared for publication should be a work of art. Your
figures should tell the entire story – the article text just fills in the blanks.
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
the showcase**.

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
