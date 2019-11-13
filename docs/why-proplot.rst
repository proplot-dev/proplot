============
Why ProPlot?
============

Efficiency
==========
Power users (e.g. people making figures for publication) generally want to change a ton of plot setting all at once. In matplotlib, this requires a bunch of one-liner setters and getters, but ProPlot has `format`. Also, matplotlib and cartopy use a bunch of verbose class names, while ProPlot helps bypass these with constructor functions like `Locator`, `Formatter`, etc. Keyword args to various funcs are implicitly passed to these constructors. I should **list the places** where these constructor functions are used.

Better subplots
===============
While matplotlib's `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself, ProPlot's `subplots` returns an `axes_grid` of axes, meant to unify these three possible return values.  `axes_grid` is a `list` subclass supporting 1D indexing (e.g. ``axs[0]``), but permits 2D indexing (e.g. ``axs[1,0]``) *just in case* the user *happened* to draw a clean 2D matrix of subplots. The `~axes_grid.__getattr__` override also means it no longer matters whether you are calling a method on an axes or a singleton `axes_grid` of axes. Finally, `axes_grid` lets `subplots` support complex arrangements of subplots -- just use 1D indexing when they don't look like a 2D matrix.

Gridspec policy
===============
While matplotlib permits arbitrarily many gridspecs per figure, ProPlot permits only *one*. When `subplots` is used, this is trivial to enforce. When `~Figure.add_subplot` is used, the figure geometry is "locked" after the first call -- although `~Figure.add_subplot` calls that divide into the existing geometry are also acceptable (for example, two square subplots above a longer rectangle subplots with the integers ``221``, ``222``, and ``212``).  This choice is not a major imposition on the user (see discussion in #50), and *considerably* simplifies gridspec adjustments, e.g. the "tight layout" adjustments. Also, the ProPlot tight layout algorithm is more accurate because GridSpecs can now have variable spacing between rows and columns (note: should add a tight layout option that keeps spaces constant even with the adjustments).

Flexible units
==============
Matplotlib uses figure-relative and axes-relative units for most things, e.g. subplot spacing. This encourages "tinkering" with meaningless numbers that change the subjective appearence when figure or axes dimensions change, because **text** and **plotted content** are specified in the physical units "points". ProPlot encourages users to always be aware of their physical units while `plot.units` allows people to use non-Imperical units and font size-relative units, like `left=2em` or `left=2Em`, which is more intuitive than figure-relative units. (note: should add `'ax'` and `'fig'` as optional units, and specify `units` as an axes method?). I should also **list** the places where ProPlot accepts physical units.

Related to the above, a common mistake I've noticed is people have their `'figure.dpi'` set too low or use a low-quality default backend, which encourages them to enlarge their figure dimensions (with `fig.set_size_inches()`) and enlarge their font sizes, but then when preparing for publication, shrink it back down! ProPlot enforces smaller font sizes + figure sizes and a higher resolution backend, so users don't need to downsize their figures when embedding in a PDF for publication -- the "points" and "inches" in the vector graphics are **correct**.

Flexible dimensions
===================
Most of the time (except when preparing publications), I don't care about the "physical size" of my figure -- just the physical size of **axes content**, e.g. how a 12 point text or 2 point thick line appears in the subplot box. This mode of thinking requires that you specify the physical dimensions of **individual sublots**, and allow figure dimensions to vary. This is not possible in matplotlib, but is in ProPlot. Likely, the reason this wasn't implemented is several matplotlib backends require figure size to be fixed -- if `fig.draw()` change the figure size, by the time the backend calls `fig.draw()` it is "too late" and you get weird bugs. I got "flexible figure dimensions" to work with the static inline backend (simple fix) and the Qt popup backend (only fixable because of a sketchy detail in its implementation; might break someday), but it is unfixable for the "notebook" and "macosx" backends, and possibly other popup backends.

Color usage tools
=================
The color usage tools are incredibly useful. Should mention how colormap names are case-insensitive, all colormaps and cycles are reversible by appending `_r` to the name. Also, "cycles" and "colormaps" are *fluid*, e.g. you can use a colormap as the color cycler for line plots. This is ProPlot's answer to the seaboarn idea of "palettes", which I think are unnecessarily convoluted -- matplotlib's "colormaps" and "property cyclers" are sufficient.

Matplotlib color handling features are also cleaned up. By default, when `extend!='both'`, matplotlib just trims the excess colors from either side of the map, but 99% of the time users want to use the full color range. ProPlot ensures the full color range is always used. Also, by default, matplotlib divides colormaps up into levels by sampling the colormap with a low resolution lookup table, which makes it hard to get boundaries just right. ProPlot samples the colormap with a high resolution lookup table and uses `BinNorm` to restrict the possible colormap indices; it makes more sense to handle **all** aspects of the value --> color conversion process, including **discretization**, to the normalizer.

Cartopy and basemap improvements
================================
Better projection handling and more sensible defaults. Part of this falls into the "use constructor functions everywhere and register verbose class objects under succinct string names" philosophy. Another is that basemap `latlon=True` should **really** be the default (99% of the time people are not working with data stored in map projection coordinates!). Same goes with cartopy `transform=ccrs.PlateCarree()` -- I think cartopy made the map projection transform the default because they wanted to copy the behavior of basemap.
