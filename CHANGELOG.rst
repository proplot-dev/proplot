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
This will be published when some major refactoring tasks are completed,
and deprecation warnings will be removed. See :pr:`89`, :pr:`109`, :pr:`110`,
and :pr:`111`.

ProPlot v0.6.0 (2020-##-##)
===========================
.. rubric:: Deprecated

- Deprecate `~proplot.axes.Axes.format` functions in favor of the
  axes-artist `~matplotlib.artist.Artist.set` override (:pr:`89`).
- Rename `width` and `height` `~proplot.subplots.subplots` keyword args to `figwidth` and `figheight` (:pr:`###`).
- Rename `aspect`, `axwidth`, and `axheight` keyword args to `refaspect`, `refwidth`, and `refheight` (:pr:`###`).
- Rename :rcraw:`subplots.pad` and :rcraw:`subplots.axpad` to
  :rcraw:`subplots.edgepad` and :rcraw:`subplots.subplotpad` (:pr:`###`).

.. rubric:: Features

- All features are now implemented with individual *setters*, like in matplotlib,
  but we still encourage using the bulk ``set`` method through documentation
  examples and by populating the ``set`` docstring (so valid arguments are no
  longer implicit).
- Users can now use `~proplot.subplots.figure` with
  `~proplot.subplots.Figure.add_subplot`
  *or* `~proplot.subplots.subplots` (:pr:`110`). This is a major improvement!
- `~proplot.subplots.GridSpec` now accepts physical units, rather than having
  `~proplot.subplots.subplots` handle the units (:pr:`110`).
- Allow "hanging" twin *x* and *y* axes as members of the `~proplot.subplots.EdgeStack`
  container. Arbitrarily many siblings are now permitted.
- Use `~proplot.subplots.GeometrySolver` for calculating various automatic layout
  stuff instead of having 1000 hidden `~proplot.subplots.Figure` methods (:pr:`110`).
- Use `~proplot.subplots.EdgeStack` class for handling
  stacks of colorbars, legends, and text (:pr:`110`).

.. rubric:: Internals

- Assignments to `~proplot.rctools.rc_configurator` are now validated, and the
  configurator is now a monkey patch of `~matplotlib.rcParams` (:pr:`109`).
- Plotting wrapper features (e.g. `~proplot.wrappers.standardize_1d`) are now
  implemented and documented on the individual methods themselves
  (e.g. `~proplot.axes.Axes.plot`; :pr:`111`).
  This is much easier for new users.
- Handle all projection keyword arguments in `~proplot.subplots.Figure.add_subplot`
  instead of `~proplot.subplots.subplots` (:pr:`110`).
- Panels, colorbars, and legends are now members of `~proplot.subplots.EdgeStack`
  stacks rather than getting inserted directly into
  the main `~proplot.subplots.GridSpec` (:pr:`110`).

ProPlot v0.5.0 (2020-02-10)
===========================
.. rubric:: Deprecated

- Remove `abcformat` from `~proplot.axes.Axes.format` (:commit:`2f295e18`).
- Rename `top` to `abovetop` in `~proplot.axes.Axes.format` (:commit:`500dd381`).
- Rename `abc.linewidth` and `title.linewidth` to ``borderwidth`` (:commit:`54eb4bee`).
- Rename `~proplot.wrappers.text_wrapper` `linewidth` and `invert` to
  `borderwidth` and `borderinvert` (:commit:`54eb4bee`).

.. rubric:: Features

- Add back `Fabio Crameri's scientific colour maps <http://www.fabiocrameri.ch/colourmaps.php>`__ (:pr:`116`).
- Permit both e.g. `locator` and `xlocator` as keyword arguments to
  `~proplot.axes.Axes.altx`, etc. (:issue:`57fab860`).
- Permit *descending* `~proplot.styletools.BinNorm` and `~proplot.styletools.LinearSegmentedNorm`
  levels (:pr:`119`).
- Permit overriding the font weight, style, and stretch in the
  `~proplot.styletools.show_fonts` table (:commit:`e8b9ee38`).
- Permit hiding "unknown" colormaps and color cycles in the
  `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  tables (:commit:`cb206f19`).

.. rubric:: Bug fixes

- Fix issue where `~proplot.styletools.show_cmaps`
  and `~proplot.styletools.show_cycles` colormap names were messed up
  (:commit:`13045599`)
- Fix issue where `~proplot.styletools.show_cmaps`
  and `~proplot.styletools.show_cycles` did not return figure instance
  (:commit:`98209e87`).
- Fix issue where user `values` passed to `~proplot.wrappers.colorbar_wrapper`
  were sometimes ignored (:commit:`fd4f8d5f`).
- Prevent formatting rightmost meridian label as ``1e-10``
  on cartopy map projections (:commit:`37fdd1eb]`).
- Support CF-time axes by fixing bug in `~proplot.wrappers.standardize_1d`
  and `~proplot.wrappers.standardize_2d` (:issue:`103`, :pr:`121`).
- Redirect to the "default" location
  when using ``legend=True`` and ``colorbar=True`` to generate on-the-fly legends
  and colorbars (:commit:`c2c5c58d`). This feature was accidentally removed.
- Let `~proplot.wrappers.colorbar_wrapper`
  accept lists of colors (:commit:`e5f11591`). This feature was accidentally removed.

.. rubric:: Internals

- Remove various unused keyword arguments (:commit:`33654a42`).
- Major improvements to the API controlling axes titles and a-b-c labels
  (:commit:`1ef7e65e`).
- Always use full names ``left``, ``right``, ``top``, and ``bottom`` instead of ``l``, ``r``,
  ``b``, and ``t``, for clarity (:commit:`1ef7e65e`).
- Improve ``GrayCycle`` colormap, is now much shorter and built from reflected
  Fabio ``GrayC`` colormaps (:commit:`5b2c7eb7`).


ProPlot v0.4.3 (2020-01-21)
===========================
.. rubric:: Deprecated

- Remove `~proplot.rctools.ipython_autoreload`,
  `~proplot.rctools.ipython_autosave`, and `~proplot.rctools.ipython_matplotlib`
  (:issue:`112`, :pr:`113`). Move inline backend configuration to a hidden
  method that gets called whenever the ``rc_configurator`` is initalized.

.. rubric:: Features

- Permit comments at the head of colormap and color files (:commit:`0ffc1d15`).
- Make `~proplot.axes.Axes.parametric` match ``plot`` autoscaling behavior
  (:commit:`ecdcba82`).

.. rubric:: Internals

- Use `~proplot.axes.Axes.colorbar` instead of `~matplotlib.axes.Axes.imshow`
  for `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  displays (:pr:`107`).

ProPlot v0.4.2 (2020-01-09)
===========================
.. rubric:: Features

- Add ``family`` keyword arg to `~proplot.styletools.show_fonts` (:pr:`106`).
- Package the `TeX Gyre <http://www.gust.org.pl/projects/e-foundry/tex-gyre>`__
  font series with ProPlot (:pr:`106`). Remove a couple other fonts.
- Put the TeX Gyre fonts at the head of the serif, sans-serif, monospace,
  cursive, and fantasy ``rcParams`` font family lists (:issue:`104`, :pr:`106`).

.. rubric:: Bug fixes

- Fix issues with Fira Math weights unrecognized by matplotlib (:pr:`106`).

ProPlot v0.4.1 (2020-01-08)
===========================
.. rubric:: Deprecation

- Change the default ``.proplotrc`` format from YAML to the ``.matplotlibrc``
  syntax (:pr:`101`).

.. rubric:: Features

- Comments (lines starting with ``#``) are now permitted in all RGB and HEX style
  colormap and cycle files (:pr:`100`).
- Break down `~proplot.styletools.show_cycles` bars into categories, just
  like `~proplot.styletools.show_cmaps` (:pr:`100`).

.. rubric:: Bug fixes

- Fix issue where `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  draw empty axes (:pr:`100`).
- Add back the :ref:`default .proplorc file <The .proplotrc file>` to docs (:pr:`101`).
  To do this, ``conf.py`` auto-generates a file in ``_static``.

.. rubric:: Internals

- Add ``geogrid.color/linewidth/etc`` and ``gridminor.color/linewidth/etc`` props
  as *children* of ``grid.color/linewidth/etc`` (:pr:`101`).
- Various `~proplot.rctools.rc_configurator` improvements, remove outdated
  global variables (:pr:`101`).
- Better error handling when loading colormap/cycle files, and calls to
  `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` now raise errors while
  calls to `~proplot.styletools.register_cmaps` and `~proplot.styletools.register_cycles`
  still issue warnings (:pr:`100`).

ProPlot v0.4.0 (2020-01-07)
===========================
.. rubric:: Deprecated

- Rename `basemap_defaults` to `~proplot.projs.basemap_kwargs` and `cartopy_projs`
  to `~proplot.projs.cartopy_names` (:commit:`431a06ce`).
- Remove ``subplots.innerspace``, ``subplots.titlespace``,
  ``subplots.xlabspace``, and ``subplots.ylabspace`` spacing arguments,
  automatically calculate default non-tight spacing using `~proplot.subplots._get_space`
  based on current tick lengths, label sizes, etc.
- Remove redundant `~proplot.rctools.use_fonts`, use ``rcParams['sans-serif']``
  precedence instead (:pr:`95`).
- `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy` no longer accept
  "scale-spec" arguments.
  Must be a function, two functions, or an axis scale instance (:pr:`96`).
- Remove `~proplot.axes.Axes` ``share[x|y]``, ``span[x|y]``, and ``align[x|y]``
  kwargs (:pr:`99`).
  These settings are now always figure-wide.
- Rename `~proplot.styletools.Cycle` ``samples`` to ``N``, rename
  `~proplot.styletools.show_colors` ``nbreak`` to ``nhues`` (:pr:`98`).

.. rubric:: Features

- Add `~proplot.styletools.LinearSegmentedColormap.from_file` static methods (:pr:`98`).
  You can now load files by passing a name to `~proplot.styletools.Colormap`.
- Add TeX Gyre Heros as open source Helvetica-alternative; this is the new default font.
  Add Fira Math as DejaVu Sans-alternative; has complete set of math characters (:pr:`95`).
- Add `xlinewidth`, `ylinewidth`, `xgridcolor`, `ygridcolor` keyword
  args to `~proplot.axes.XYAxes.format` (:pr:`95`).
- Add getters and setters for various `~proplot.subplots.Figure` settings
  like ``share[x|y]``, ``span[x|y]``, and ``align[x|y]`` (:pr:`99`).
- Let `~proplot.axes.Axes.twinx`, `~proplot.axes.Axes.twiny`,
  `~proplot.axes.Axes.altx`, and `~proplot.axes.Axes.alty` accept
  `~proplot.axes.XYAxes.format` keyword args just like
  `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy` (:pr:`99`).
- Add `~proplot.subplots.Figure` ``fallback_to_cm`` kwarg. This is used by
  `~proplot.styletools.show_fonts` to show dummy glyphs to clearly illustrate when fonts are
  missing characters, but preserve graceful fallback for end user.
- Improve `~proplot.projs.Proj` constructor function. It now accepts
  `~cartopy.crs.Projection` and `~mpl_toolkits.basemap.Basemap` instances, just like other
  constructor functions, and returns only the projection instance (:pr:`92`).
- `~proplot.rctools.rc` `~proplot.rctools.rc_configurator.__getitem__` always
  returns the setting. To get context block-restricted settings, you must explicitly pass
  ``context=True`` to `~proplot.rctools.rc_configurator.get`, `~proplot.rctools.rc_configurator.fill`,
  or `~proplot.rctools.rc_configurator.category` (:pr:`91`).

.. rubric:: Bug fixes

- Fix `~proplot.rctools.rc_configurator.context` bug (:issue:`80` and :pr:`91`).
- Fix issues with `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy`
  with non-linear parent scales (:pr:`96`).
- Ignore TTC fonts because they cannot be saved in EPS/PDF figures (:issue:`94` and :pr:`95`).
- Do not try to use Helvetica Neue because "thin" font style is read as regular (:issue:`94` and :pr:`95`).

.. rubric:: Documentation

- Use the imperative mood for docstring summaries (:pr:`92`).
- Fix `~proplot.styletools.show_cycles` bug (:pr:`90`) and show cycles using colorbars
  rather than lines (:pr:`98`).

.. rubric:: Internals

- Define `~proplot.rctools.rc` default values with inline dictionaries rather than
  with a default ``.proplotrc`` file, change the auto-generated user ``.proplotrc``
  (:pr:`91`).
- Remove useless `panel_kw` keyword arg from `~proplot.wrappers.legend_wrapper` and
  `~proplot.wrappers.colorbar_wrapper` (:pr:`91`). Remove `wflush`, `hflush`,
  and `flush` keyword args from `~proplot.subplots.subplots` that should have been
  removed long ago.

ProPlot v0.3.1 (2019-12-16)
===========================
.. rubric:: Bug fixes

- Fix issue where custom fonts were not synced (:commit:`a1b47b4c`).
- Fix issue with latest versions of matplotlib where ``%matplotlib inline``
  fails *silently* so the backend is not instantiated (:commit:`cc39dc56`).

ProPlot v0.3.0 (2019-12-15)
===========================
.. rubric:: Deprecated

- Remove ``'Moisture'`` colormap (:commit:`cf8952b1`).

.. rubric:: Features

- Add `~proplot.styletools.use_font`, only sync Google Fonts fonts (:pr:`87`).
- New ``'DryWet'`` colormap is colorblind friendly (:commit:`0280e266`).
- Permit shifting arbitrary colormaps by ``180`` degrees by appending the
  name with ``'_shifted'``, just like ``'_r'`` (:commit:`e2e2b2c7`).

.. rubric:: Bug fixes

- Add brute force workaround for saving colormaps with
  *callable* segmentdata (:commit:`8201a806`).
- Fix issue with latest versions of matplotlib where ``%matplotlib inline``
  fails *silently* so the backend is not instantiated (:commit:`cc39dc56`).
- Fix `~proplot.styletools.LinearSegmentedColormap.shifted` when `shift` is
  not ``180`` (:commit:`e2e2b2c7`).
- Save the ``cyclic`` and ``gamma`` attributes in JSON files too (:commit:`8201a806`).

.. rubric:: Documentation

- Cleanup notebooks, especially the colormaps demo (e.g. :commit:`952d4cb3`).

.. rubric:: Internals

- Change `~time.clock` to `~time.perf_counter` (:pr:`86`).

ProPlot v0.2.7 (2019-12-09)
===========================

.. rubric:: Bug fixes

- Fix issue where `~proplot.styletools.AutoFormatter` logarithmic scale
  points are incorrect (:commit:`9b164733`).

.. rubric:: Documentation

- Improve :ref:`Configuring proplot` documentation (:commit:`9d50719b`).

.. rubric:: Internals

- Remove `prefix`, `suffix`, and `negpos` keyword args from
  `~proplot.styletools.SimpleFormatter`, remove `precision` keyword arg from
  `~proplot.styletools.AutoFormatter` (:commit:`8520e363`).
- Make ``'deglat'``, ``'deglon'``, ``'lat'``, ``'lon'``, and ``'deg'`` instances
  of `~proplot.styletools.AutoFormatter` instead of `~proplot.styletools.SimpleFormatter`
  (:commit:`8520e363`). The latter should just be used for contours.

ProPlot v0.2.6 (2019-12-08)
===========================
.. rubric:: Bug fixes

- Fix issue where twin axes are drawn *twice* (:commit:`56145122`).


ProPlot v0.2.5 (2019-12-07)
===========================
.. rubric:: Features

- Much better `~proplot.axistools.CutoffScale` algorithm, permit arbitrary
  cutoffs (:pr:`83`).

ProPlot v0.2.4 (2019-12-07)
===========================
.. rubric:: Deprecated

- Rename `ColorCacheDict` to `~proplot.styletools.ColorDict` (:commit:`aee7d1be`).
- Rename `colors` to `~proplot.styletools.Colors` (:commit:`aee7d1be`)
- Remove `fonts_system` and `fonts_proplot`, rename `colordict` to
  `~proplot.styletools.colors`, make top-level variables
  more robust (:commit:`861583f8`).

.. rubric:: Documentation

- Params table for `~proplot.styletools.show_fonts` (:commit:`861583f8`).

.. rubric:: Internals

- Improvements to `~proplot.styletools.register_colors`.

ProPlot v0.2.3 (2019-12-05)
===========================
.. rubric:: Bug fixes

- Fix issue with overlapping gridlines using monkey patches on gridliner
  instances (:commit:`8960ebdc`).
- Fix issue where auto colorbar labels are not applied when
  ``globe=True`` (:commit:`ecb3c899`).
- More sensible zorder for gridlines (:commit:`90d94e55`).
- Fix issue where customized super title settings are overridden when
  new axes are created (:commit:`35cb21f2`).

.. rubric:: Documentation

- Organize ipython notebook documentation (:commit:`35cb21f2`).

.. rubric:: Internals

- Major cleanup of the `~proplot.wrappers.colorbar_wrapper` source code, handle
  minor ticks using the builtin matplotlib API just like major ticks (:commit:`b9976220`).

ProPlot v0.2.2 (2019-12-04)
===========================
.. rubric:: Deprecated

- Rename `~proplot.subplots.axes_grid` to `~proplot.subplots.subplot_grid` (:commit:`ac14e9dd`).

.. rubric:: Bug fixes

- Fix shared *x* and *y* axis bugs (:commit:`ac14e9dd`).

.. rubric:: Documentation

- Make notebook examples PEP8 compliant (:commit:`97f5ffd4`). Much more readable now.

ProPlot v0.2.1 (2019-12-02)
===========================
.. rubric:: Deprecated

- Rename `autoreload_setup`, `autosave_setup`, and `matplotlib_setup` to
  `~proplot.rctools.ipython_autoreload`, `~proplot.rctools.ipython_autosave`, and `~proplot.rctools.ipython_matplotlib`, respectively (:commit:`84e80c1e`).

ProPlot v0.2.0 (2019-12-02)
===========================
.. rubric:: Deprecated

- Remove the ``nbsetup`` rc setting in favor of separate ``autosave``, ``autoreload``,
  and ``matplotlib`` settings for triggering the respective ``%`` magic commands.
  (:commit:`3a622887`; ``nbsetup`` is still accepted but no longer documented).
- Rename the ``format`` rc setting in favor of the ``inlinefmt`` setting
  (:commit:`3a622887`; ``format`` is still accepted but no longer documented).
- Rename ``FlexibleGridSpec`` and ``FlexibleSubplotSpec`` to ``GridSpec``
  and ``SubplotSpec`` (:commit:`3a622887`; until :pr:`110` is merged it is impossible
  to use these manually, so this won't bother anyone).

.. rubric:: Features

- Support manual resizing for all backends, including ``osx`` and ``qt`` (:commit:`3a622887`).

.. rubric:: Bug fixes

- Disable automatic resizing for the ``nbAgg`` interactive inline backend. Found no
  suitable workaround (:commit:`3a622887`).

.. rubric:: Internals

- Organize the ``rc`` documentation and the default ``.proplotrc`` file (:commit:`3a622887`).
- Rename ``rcParamsCustom`` to ``rcParamsLong``
  (:commit:`3a622887`; this is inaccessible to the user).
- When calling ``fig.canvas.print_figure()`` on a stale figure, call ``fig.canvas.draw()``
  first. May be overkill for `~matplotlib.figure.Figure.savefig` but critical for
  correctly displaying already-drawn notebook figures.

ProPlot v0.1.0 (2019-12-01)
===========================
.. rubric:: Internals

- Include `flake8` in Travis CI testing (:commit:`8743b857`).
- Enforce source code PEP8 compliance (:commit:`78da51a7`).
- Use pre-commit for all future commits (:commit:`e14f6809`).
- Implement tight layout stuff with canvas monkey patches (:commit:`67221d10`).
  ProPlot now works for arbitrary backends, not just inline and qt.

.. rubric:: Documentation

- Various `RTD bugfixes <https://github.com/readthedocs/readthedocs.org/issues/6412>`__ (e.g. :commit:`37633a4c`).

ProPlot v0.0.0 (2019-11-27)
===========================

The first version released on `PyPi <https://pypi.org/project/proplot/>`__.

.. _`Luke Davis`: https://github.com/lukelbd
.. _`Riley X. Brady`: https://github.com/bradyrx
