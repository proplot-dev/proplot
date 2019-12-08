..
  Valid subsections:
  - Deprecated
  - Features
  - Bug fixes
  - Internals
  - Documentation

=================
Changelog history
=================

ProPlot v1.0.0 (2020-##-##)
===========================
This will be published when some major refactoring tasks are completed.
See :pr:`45`, :pr:`46`, and :pr:`50`.

ProPlot v0.3.0 (2020-01-XX)
===========================
Features
--------
- Users can now use `~proplot.subplots.figure` with `~proplot.subplots.Figure.add_subplot`
  *or* `~proplot.subplots.subplots` (:pr:`50`). This is a major improvement!
- `~proplot.subplots.GridSpec` now accepts physical units, rather than having
  `~proplot.subplots.subplots` handle the units (:pr:`50`).
- Add `xlinewidth`, `ylinewidth`, `xgridcolor`, `ygridcolor` keyword
  args to `~proplot.axes.XYAxes.format` (:pr:`50`).
- Allow "hanging" twin *x* and *y* axes as members of the `~proplot.subplots.EdgeStack`
  container. Arbitrarily many siblings are now permitted.
- Use `~proplot.subplots.GeometrySolver` for calculating various automatic layout
  stuff instead of having 1000 hidden `~proplot.subplots.Figure` methods (:pr:`50`).
- Use `~proplot.subplots.EdgeStack` class for handling
  stacks of colorbars, legends, and text (:pr:`50`).

Bug fixes
---------
- Fix `~proplot.rctools.rc_configurator.context` fatal bug (:issue:`80`).

Deprecated
----------
- Remove ``subplots.innerspace``, ``subplots.titlespace``,
  ``subplots.xlabspace``, and ``subplots.ylabspace`` spacing arguments,
  automatically calculate default non-tight spacing using `~proplot.subplots._get_space`
  based on current tick lengths, label sizes, etc.

Internals
---------
- Handle all projection keyword arguments in `~proplot.subplots.Figure.add_subplot`
  instead of `~proplot.subplots.subplots` (:pr:`50`).
- Panels, colorbars, and legends are now members of `~proplot.subplots.EdgeStack`
  stacks rather than getting inserted directly into
  the main `~proplot.subplots.GridSpec` (:pr:`50`).
- Define `~proplot.rctools.rc` default values with inline dictionaries rather than
  with a default ``.proplotrc`` file, change the auto-generated user ``.proplotrc``
  (:pr:`50`).

ProPlot v0.2.6 (2019-12-08)
===========================
Bug fixes
---------
- Fix issue where twin axes are drawn *twice* (:commit:`56145122`).


ProPlot v0.2.5 (2019-12-07)
===========================
Features
--------
- Much better `~proplot.axistools.CutoffScale` algorithm, permit arbitrary
  cutoffs (:pr:`83`).

ProPlot v0.2.4 (2019-12-07)
===========================
Deprecated
----------
- Rename `ColorCacheDict` to `~proplot.styletools.ColorDict` (:commit:`aee7d1be`).
- Rename `colors` to `~proplot.styletools.Colors` (:commit:`aee7d1be`)
- Remove `fonts_system` and `fonts_proplot`, rename `colordict` to
  `~proplot.styletools.colors`, make top-level variables
  more robust (:commit:`861583f8`).

Documentation
-------------
- Params table for `~proplot.styletools.show_fonts` (:commit:`861583f8`).

Internals
---------
- Improvements to `~proplot.styletools.register_colors`.

ProPlot v0.2.3 (2019-12-05)
===========================
Bug fixes
---------
- Fix issue with overlapping gridlines (:commit:`8960ebdc`).
- Fix issue where auto colorbar labels are not applied when ``globe=True`` (:commit:`ecb3c899`).
- More sensible zorder for gridlines (:commit:`90d94e55`).
- Fix issue where customized super title settings are overridden when
  new axes are created (:commit:`35cb21f2`).

Documentation
-------------
- Organize ipython notebook documentation (:commit:`35cb21f2`).

Internals
---------
- Major cleanup of the `~proplot.wrappers.colorbar_wrapper` source code, handle
  minor ticks using the builtin matplotlib API just like major ticks (:commit:`b9976220`).

ProPlot v0.2.2 (2019-12-04)
===========================
Bug fixes
---------
- Fix shared *x* and *y* axis bugs (:commit:`ac14e9dd`).

Deprecated
----------
- Rename `~proplot.subplots.axes_grid` to `~proplot.subplots.subplot_grid` (:commit:`ac14e9dd`).

Documentation
-------------
- Make notebook examples PEP8 compliant (:commit:`97f5ffd4`). Much more readable now.

ProPlot v0.2.1 (2019-12-02)
===========================
Deprecated
----------
- Rename `autoreload_setup`, `autosave_setup`, and `matplotlib_setup` to
  `~proplot.rctools.ipython_autoreload`, `~proplot.rctools.ipython_autosave`, and `~proplot.rctools.ipython_matplotlib`, respectively (:commit:`84e80c1e`).

ProPlot v0.2.0 (2019-12-02)
===========================
Features
--------
- Support manual resizing for all backends, including ``osx`` and ``qt`` (:commit:`3a622887`).

Bug fixes
---------
- Disable automatic resizing for the ``nbAgg`` interactive inline backend. Found no
  suitable workaround (:commit:`3a622887`).

Deprecated
----------
- Remove the ``nbsetup`` rc setting in favor of separate ``autosave``, ``autoreload``,
  and ``matplotlib`` settings for triggering the respective ``%`` magic commands.
  (:commit:`3a622887`; ``nbsetup`` is still accepted but no longer documented).
- Rename the ``format`` rc setting in favor of the ``inlinefmt`` setting
  (:commit:`3a622887`; ``format`` is still accepted but no longer documented).
- Rename ``FlexibleGridSpec`` and ``FlexibleSubplotSpec`` to ``GridSpec``
  and ``SubplotSpec`` (:commit:`3a622887`; until :pr:`50` is merged it is impossible
  to use these manually, so this won't bother anyone).

Internals
---------
- Organize the ``rc`` documentation and the default ``.proplotrc`` file (:commit:`3a622887`).
- Rename ``rcParamsCustom`` to ``rcParamsLong``
  (:commit:`3a622887`; this is inaccessible to the user).
- When calling ``fig.canvas.print_figure()`` on a stale figure, call ``fig.canvas.draw()``
  first. May be overkill for `~matplotlib.figure.Figure.savefig` but critical for
  correctly displaying already-drawn notebook figures.

ProPlot v0.1.0 (2019-12-01)
===========================
Internals
---------
- Include `flake8` in Travis CI testing (:commit:`8743b857`).
- Enforce source code PEP8 compliance (:commit:`78da51a7`).
- Use pre-commit for all future commits (:commit:`e14f6809`).
- Implement tight layout stuff with canvas monkey patches (:commit:`67221d10`).
  ProPlot now works for arbitrary backends, not just inline and qt.

Documentation
-------------
- Various `RTD bugfixes <https://github.com/readthedocs/readthedocs.org/issues/6412>`__ (e.g. :commit:`37633a4c`).

ProPlot v0.0.0 (2019-11-27)
===========================

The first version released on `PyPi <https://pypi.org/project/proplot/>`__.

.. _`Luke Davis`: https://github.com/lukelbd
.. _`Riley X. Brady`: https://github.com/bradyrx
