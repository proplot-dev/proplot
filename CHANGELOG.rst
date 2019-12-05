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

ProPlot v0.2.X (2019-12-02)
===========================
Features
--------
- Support manual resizing for all backends, including ``osx`` and ``qt`` (:commit:`3a622887`).

Bug fixes
---------
- Fix issue with overlapping gridlines (:commit:`8960ebdc`).
- Fix issue where auto colorbar labels are not applied when ``globe=True`` (:commit:`ecb3c899`)
- More sensible zorder for gridlines (:commit:`90d94e55`).
- Fix shared *x* and *y* axis bugs (:commit:`ac14e9dd`).
- Disable automatic resizing for the ``nbAgg`` interactive inline backend. Found no
  suitable workaround (:commit:`3a622887`).
- Fix issue where customized super title settings are overridden when
  new axes are created (:commit:`35cb21f2`)

Deprecated
----------
- Rename `~proplot.subplots.axes_grid` to `~proplot.subplots.subplot_grid` (:commit:`ac14e9dd`).
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
- Rename ``rcParamsCustom`` to ``rcParamsLong`` (this is inaccessible to the user).
- Organize the ``rc`` documentation and the default ``.proplotrc`` file.
- When calling ``fig.canvas.print_figure()`` on a stale figure, call ``fig.canvas.draw()``
  first. May be overkill for `~matplotlib.figure.Figure.savefig` but critical for
  correctly displaying already-drawn notebook figures.

Documentation
-------------
- Make notebook examples PEP8 compliant (:commit:`97f5ffd4`). Much more readable now.
- Clean up documentation (:commit:`35cb21f2`).

ProPlot v0.1.X (2019-12-01)
===========================
Internals
---------
- Include `flake8` in Travis CI testing (:commit:`8743b857`).
- Enforce source code PEP8 compliance (:commit:`78da51a7`).
- Use pre-commit for all future commits (:commit:`e14f6809`).
- Implement tight layout stuff with canvas monkey patches (:commit:`67221d10`).
  This is more robust to different backends.

Documentation
-------------
- Various `RTD bugfixes <https://github.com/readthedocs/readthedocs.org/issues/6412>`__ (e.g. :commit:`37633a4c`).

ProPlot v0.0.0 (2019-11-27)
===========================

The first version released on `PyPi <https://pypi.org/project/proplot/>`__.

.. _`Luke Davis`: https://github.com/lukelbd
.. _`Riley X. Brady`: https://github.com/bradyrx
