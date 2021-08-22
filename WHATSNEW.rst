..
  Valid rubrics:
  - Deprecated
  - Features
  - Bug fixes
  - Internals
  - Documentation

.. _whats_new:

===========
What's new?
===========

This page lists the API changes with each version. Please note that
when functions and keywords are renamed, they are not removed -- requesting
the old name will issue a warning rather than an error. Authors are shown next to
each change. Where not indicated, `Luke Davis`_ was the author.

See the :ref:`author page <authors>` for a list of contributors, and see
the :ref:`contribution guide <contributions>` if you are interested in
submitting your own changes.

ProPlot v1.0.0 (2022-##-##)
===========================

This will be published when more comprehensive testing is completed
and stability is improved.

ProPlot v0.8.1 (2021-08-22)
===========================

Features
--------

* Add `~proplot.colors.PerceptualColormap.from_list` ``adjust_grays`` option
  (enabled by default) to help make diverging colormaps with an intermediate
  hueless white, gray, or black color (:commit:`2e8cb495`).
* Add the axis sharing level ``4`` or ``'all'`` to share the limits, scales,
  and tick labels between axes not in the same row/column (:commit:`73f355a2`).
* Allow adding contours to `legend` by interpreting `label` keyword and using
  central handle from ``ContourSet.legend_elements`` (:commit:`26bc77a4`).
* Extend mixed auto-manual legend label input (e.g. ``labels=[None, 'override']``)
  to case where legend handles are automatically retrieved from the axes
  rather than manually passed to ``legend()`` (:commit:`26bc77a4`).
* Add `inlinelabels` option to `~proplot.axes.GeoAxes.format` to set both
  ``loninline=True`` and ``latinline=True`` at once, and change the
  :rcraw:`grid.loninline` and :rcraw:`grid.latinline` settings to the
  single :rcraw:`grid.inlinelabels` (consistent with :rcraw:`grid.rotatelabels`
  and :rcraw:`grid.dmslabels`) (:commit:`560ed978`).

Bug fixes
---------

* Fix regression where ``np.std`` and ``np.percentile`` no longer
  ignore NaN values (:issue:`257`, :commit:`d1906fce`).
* Fix regression where ``legend()`` cannot be called without
  the input labels (:issue:`188`, :commit:`fdd53a6c`).
* Fix issue where edge colors of area plots with ``negpos=True``
  cannot be changed (:commit:`bb50dea4`).
* Fix issue where `legend` `order` keyword arg is ignored and default is
  changed back to ``'F'`` (:commit:`###`).
* Fix issues where ``setup_matplotlib`` is not called for pint quantity
  input and column iteration of 2D input to 1D funcs fails (:commit:`e57d238e`).
* Fix issue where pint quantity *x* and *y* coordinates fail when passing
  as pcolor centers or when :rcraw:`cmap.inbounds` enabled (:commit:`fd76af3a`).
* Fix issue where pint quantity *z* data do not have units stripped
  unless in xarray dataarray (:commit:`aadc65f9`).
* Fix issue where making single-color contour plots creates just one contour by
  making default ``levels`` count independent from `colors` (:commit:`63eaf10e`).
* Fix issue where common legend handle properties cannot be overridden due to
  searching for ``collection`` props rather than ``line`` props (:commit:`26bc77a4`).
* Fix issue where title/abc padding is overwritten in the presence of top panels
  and make title deflection to top panels generally more robust (:commit:`d27d05cf`).
* Fix issues with the ``%qt`` backend using ``forward=False``
  during subplot additions (:issue:`244`, :commit:`ac12bbc2`)
* Fix issue where ``%matpolotlib notebook`` and ``%matplotlib widget`` display
  unusable/cutoff figure previews by fixing the figure size at creation time and
  issuing one-time warning if size was not fixed explicitly (:commit:`88fc2868`).

Documentation
-------------

* Make docstring utils explicitly private and convert `_snippets` dictionary to
  callable dictionary-like `_SnippetsManager` instance (:commit:`b73fe9e3`). This
  helps prevent bug where assigned snippets have unfilled ``%(snippet)s`` markers.

ProPlot v0.8.0 (2021-08-18)
===========================

Deprecated
----------

* Numbers passed to `pad`, `wpad`, `hpad`, `space`, `wspace`, `hspace`, `left`, `right`,
  `top`, and `bottom` are now interpreted as em-widths instead of inches
  (:commit:`20502345`). Unfortunately this is a major breaking change that cannot be
  "gently" phased in with warnings, but this will be much more convenient going forward.
* Make default reference subplot size, panel widths, colorbar widths independent of
  :rcraw:`font.size` (:commit:`a50d5264`). Default space size should definitely sync
  with font size, since larger fonts produce larger labels between subplots, but the
  same reasoning does not apply for subplot size.
* Interpret ``sharex/sharey=True`` as ``3`` (i.e., "turn all sharing on") instead
  of ``1`` (integer conversion of ``True``) (:issue:`51967ce3`). This is more
  intuitive and matches convention elsewhere. Also allow specifying level 1 with
  ``'labels'`` and level 2 with ``'limits'``.
* Rename `~proplot.ui.SubplotsContainer` to simpler `~proplot.figure.SubplotGrid`
  and move definition to ``figure.py`` (:commit:`51967ce3`).
* Deprecate arbitrary ``__getattr__`` override for `~proplot.figure.SubplotGrid`
  (:commit:`51967ce3`). Instead have dedicated ``format``, ``colorbar``, ``legend``,
  ``[alt|dual|twin][xy]``, ``panel[_axes]``, and ``inset[_axes]`` methods.
* Deprecate scattershot `~proplot.figure.Figure` immutable/documented
  properties (:commit:`51967ce3`). Really not necessary.
* Remove ancient deprecated getters and setters for ``sharex``, ``spanx``, etc.
  once used with figure objects (:commit:`51967ce3`).
* Make `~proplot.axes.CartopyAxes` and `~proplot.axes.BasemapAxes` private and
  remove documentation (:commit:`25e759b0`). These classes are just for internal
  implementation of different cartographic "backends" -- behavior of public
  methods is the same for both. Instead just document `proplot.axes.GeoAxes`.
* Rename `~proplot.axes.Axes3D` to `~proplot.axes.ThreeAxes` so that class name
  fits more nicely amongst other class names (:commit:`30a112bd`).
* Sync 3D axes figure background color with axes background to avoid weird
  misaligned white square behind axes (:commit:`30a112bd`).
* Rename `~proplot.config.RcConfigurator` to `~proplot.config.Configurator`
  (:commit:`5626bc88`). Previous name was redundant and needlessly verbose
  (the ``c`` in ``rc`` already stands for "configuration"...).
* Add :rcraw:`leftlabel.rotation`, :rcraw:`toplabel.rotation`,
  :rcraw:`rightlabel.rotation`, :rcraw:`bottomlabel.rotation` settings, and make
  default row label rotation match y label rotation (:commit:`bae85113`).
* Rename ``image`` category settings to :rcraw:`cmap.inbounds`,
  :rcraw:`cmap.discrete`, :rcraw:`cmap.edgefix`, :rcraw:`cmap.levels`, and
  :rcraw:`cmap.lut` (:commit:`a50d5264`).
* Rename confusing :rcraw:`text.labelsize` and :rcraw:`text.titlesize` settings
  to clearer :rcraw:`font.smallsize` and :rcraw:`font.largesize` with shorthands
  :rcraw:`font.small` and :rcraw:`font.large` (analogous to :rcraw:`font.size`)
  (:pr:`a50d5264`). Previous names were bad because "label size" applies to more than
  just axis or tick labels and "title size" applies to more than just axes titles.
* Rename :rcraw:`tick.ratio` to :rcraw:`tick.widthratio` and add missing
  :rcraw:`tick.width` setting (:commit:`a50d5264`).
* Rename vague shorthands :rcraw:`alpha` and :rcraw:`facecolor` back to native
  :rcraw:`axes.alpha` and :rcraw:`axes.facecolor` and rename :rcraw:`linewidth`
  and :rcraw:`color` to :rcraw:`meta.width` and :rcraw:`meta.color`
  (:commit:`41b5e400`). However axes can still be updated by passing `alpha`,
  `linewidth`, `facecolor`, and `edgecolor` to ``format``, and now ``format`` supports
  *arbitrary* patch artist settings and aliases like `lw`, `ec`, `fc`, `hatch`, etc.
* Treat :rcraw:`tick.label` and :rcraw:`grid.label` font size, color, and weight
  settings as *synonyms* (:commit:`a50d5264`). In general the tick vs. grid distinction
  is not meaningful for text labels. However we often want different padding so still
  allow :rcraw:`tick.labelpad` and :rcraw:`grid.labelpad` to be distinct.
* Rename `proplot.config.Configurator.load_file` to `proplot.config.Configurator.load`
  in order to match ``save`` (:commit:`1769d349`).
* Change the default `~proplot.config.Configurator` save location from the home
  directory to the *current directory* and change the default filename to
  ``proplotrc`` (without the leading dot) (:commit:`41b5e400`).
* Rename obscure `LinearSegmentedColormap`, `PerceptuallyUniformColormap`, and
  `ListedColormap` to more intuitive/succinct `~proplot.colors.ContinuousColormap`,
  `~proplot.colors.PerceptualColormap`, and `~proplot.colors.DiscreteColormap`
  (:commit:`ade787f9`). Important due to the "qualitative colormap" behaviors triggered
  when a `~proplot.colors.DiscreteColormap` is passed to plot commands (see features).
* Following above change, rename `LinearSegmentedNorm` to simpler `SegmentedNorm`,
  rename `~proplot.constructor.Colormap` argument `to_listed` to `discrete`,
  change `listmode` options from ``'listed'``, ``'linear'`` to ``'discrete'``,
  ``'continuous'``, and add `filemode` option (:commit:`ade787f9`, :commit:`5ccd6c01`).
* Allow omitting the colormap name when instantiating colormap classes or using
  class methods like ``from_list`` (:commit:`ade787f9`). This is more intuitive.
* Capture `colors` passed to commands like ``contour`` and ``pcolor`` and use
  it to build qualitative `~proplot.colors.DiscreteColormap` maps (:commit:`6382cf91`).
  This matches the behavior of xarray plotting utilities. No longer use `color`
  to change "edge color" of filled contours/grid boxes.
* Set default linewidth to 0.25 when adding "edges" to filled contours
  (:commit:`6382cf91`). This matches matplotlib behavior when passing
  edgecolor to a ``pcolor`` command.
* Deprecate ``boxes`` and ``violins`` shorthands in favor of singular
  `~proplot.axes.PlotAxes.box` and `~proplot.axes.PlotAxes.violin`
  (:commit:`6382cf91`). This feel analogous to existing ``bar`` and ``barh``.
* Treat 2D ``scatter`` arguments by iterating over columns and default-styling each
  column with the property cycle rather than unraveling 2D arguments into 1D
  arrays (:commit:`6382cf91`). Can also iterate over ``s`` and ``c`` columns.
* Rename the confusingly-capitalized `~proplot.constructor.Colors` to
  `~proplot.utils.get_colors` and move to ``utils.py`` (:commit:`51d480da`). This
  is not a "class constructor" -- it just returns lists of colors.
* Rename the ``show`` function keyword `categories` to `include`,
  consistent with the new `ignore` keyword (:commit:`c45d5fa1`).
* Change default :rcraw:`legend.facecolor` to white instead of inheriting from
  axes background (:commit:`6382cf91`). Also set default :rcraw:`legend.edgecolor`
  to :rcraw:`meta.color` (black by default) and have `legend` read from rc
  settings rather than setting default `legend` input arguments.
* Control edge width for legend frames with `ew` or `edgewidth` rather than
  `lw` and `linewidth` to avoid conflict with feature that permits modifying
  legend handle properties (:commit:`6382cf91`).

Features
--------

* Dynamically add classes that are "registered" by contructor functions
  to the top-level namespace (:commit:`4382a1b1`). This is consistent with behavior
  of importing custom-ProPlot tickers, norms, etc. to top-level namespace.
  Now e.g. ``pplt.MultipleLocator`` or ``pplt.LogNorm`` are allowed.
* Allow creating subplots with `~proplot.ui.figure` and either (1) subsequently
  calling `~proplot.figure.Fiugure.subplots` or (2) passing integers or subplot specs
  generated by `~proplot.gridspec.GridSpec` to `~proplot.figure.Figure.add_subplot`
  (:commit:`51967ce3`). This is convenient for complex grids or mixed proj types.
* Add consistent/intuitive aliases `~proplot.figure.Figure.subplot` and
  `~proplot.figure.Figure.add_subplots` for native matplotlib commands
  `~proplot.figure.Figure.add_subplot` and `~proplot.figure.Figure.subplots`
  (:commit:`51967ce3`).
* Add `~proplot.figure.Figure.subplotgrid` property to access a
  `~proplot.figure.SubplotGrid` after drawing subplots one-by-one
  (:commit:`fb83384f`).
* Implement physical-units `left`, `right`, `top`, `bottom`, `wspace`, and `hspace`
  spaces directly on the `~proplot.gridspec.GridSpec` rather than externally
  (:commit:`20502345`). Now absolute spaces are always preserved when figure size
  changes even if tight layout is disabled.
* Have `~proplot.gridspec.GridSpec` directly handle "panel slots" (:commit:`20502345`).
  Adding panels to a figure adds row or column "panel slots" to the gridspec and
  subsequently indexing the gridspec ignores those slots.
* Add tight layout "padding" arguments to `~proplot.gridspec.GridSpec` and add gridspec
  parameters as optional arguments to `~proplot.figure.Figure` (:commit:`20502345`).
  When a gridspec is added to the figure the arguments are passed to the gridspec. This
  replaces matplotlib's `subplotpars` and ``subplots_adjust``.
* Allow variable tight layout padding between subplot panels using `wpad` and
  `hpad`, analogous to `wspace` and `hspace` (:commit:`20502345`). Previously
  this was fixed at :rcraw:`subplots.innerpad`.
* Add `pad` keyword to `legend`, `colorbar`, and `panel` that controls local
  tight layout padding, analogous to `space` (:commit:`20502345`). Previously this
  was fixed at :rcraw:`subplots.panelpad`.
* Ensure `wequal` and `hequal` only apply to the main subplot rows and columns;
  always ignore panel and colorbar spaces (:commit:`20502345`).
* Improve default behavior in presence of 'outer' colorbars + legends when
  :rcraw:`subplots.tight` is disabled (:commit:`20502345`).
* Add a `~proplot.figure.Figure.format` method for formatting every subplot in
  the figure when you don't have a ``SubplotGrid`` available (:commit:`20502345`).
  Also move internal implementation of figure-wide settings there. Figure-wide
  settings like `suptitle` can still be updated from ``Axes.format``.
* Permit mutability of `~proplot.figure.SubplotGrid` (:commit:`51967ce3`).
  Power users may want to manipulate their own grids.
* Permit 2d indexing of `~proplot.figure.SubplotGrid` with arbitrary gridspec
  geometry by looking up subplotspec indices (:commit:`51967ce3`). Previously 2d
  indexing of ``SubplotGrid`` with complex geometry would just return a wrong result.
* Issue warning message when users try ``fig.subplots_adjust()`` or
  ``pplt.figure(subplotpars=SubplotParams)`` and auto-disable and warn when
  matplotlib "tight layout" rc settings are toggled (:commit:`51967ce3`).
* Add nicer string representations of figures, gridspecs, subplotspecs, and
  axes clearly showing the geometry and layout (:commit:`51967ce3`, :commit:`6382cf91`).
* Rename setting :rcraw:`abc.style` to :rcraw:`abc` (:commit:`a50d5264`). Setting this
  to ``False`` still "turns off" labels, setting to ``True`` "turns on" labels with
  the default style ``'a'``, and setting to a string "turns on" labels with this style.
* Set default location for new axes panels to ``'right'``, allowing for empty
  ``ax.panel_axes()`` calls (:commit:`51967ce3`).
* Convert valid keyword arguments to positional arguments for virtually all
  plotting functions rather than a subset (:commit:`6382cf91`). This expands the
  use of the `data` keyword and permits a seaborn-like workflow (for example,
  ``ax.plot(x='x_key', y='y_key', data=xarray_dataset)``).
* Support `pint.Quantity` arguments by auto-applying ``setup_matplotlib`` with
  the quantity's unit registry when a quantity is passed (:commit:`6382cf91`).
* Support `pint.Quantity` input for *z* coordinates (e.g., to ``ax.contourf``)
  by stripping the units to prevent warning (:commit:`6382cf91`).
* Support `xarray.DataArray` arguments containing `pint.Quantity` arrays by
  accessing ``data`` rather than accessing ``.values`` (:commit:`6382cf91`).
* Apply `pint.Quantity` default unit labels to plots by formatting the units
  with the new :rcraw:`unitformat` setting (:commit:`6382cf91`).
* Add :rc:`cmap.sequential`, :rc:`cmap.diverging`, :rc:`cmap.cyclic`, and
  :rc:`cmap.qualitative` settings to control the default sequential, diverging,
  cyclic, and qualitative cmaps, and add boolean `sequential`, `diverging`, `cyclic`,
  and `qualitative` keywords to select corresponding default cmaps (:commit:`6382cf91`).
* Add special qualitative cmap handling when ``colors=colors``, ``qualitative=True``,
  or ``cmap=pcolors.DiscreteColormap(...)`` -- always apply ``DiscreteNorm`` (ignore
  and warn if user passed ``discrete=False``), truncate or wrap colors if there are too
  many/not enough for the levels, and add default extremes with ``set_under`` or
  ``set_over`` depending on user `extend` (:commit:`6382cf91`).
* Select :rcraw:`cmap.diverging` and apply `~proplot.colors.DivergingNorm` automatically
  based on input data, similar to xarray and seaborn (:commit:`6382cf91`). This is
  controlled with `autodiverging` and the :rcraw:`cmap.autodiverging` setting. It is
  also disabled when a cmap is explicitly passed (unless it is a known diverging cmap).
* Add `robust` keyword argument and :rc:`cmap.robust` setting to ignore
  outliers when selecting auto colormap ranges (:issue:`6382cf91`). It can take the
  value ``True``, a percentile range, or a 2-tuple percentile interval.
* Improve matplotlib-proplot colormap translation by converting
  `matplotlib.colors.ListedColormap` to `proplot.colors.DiscreteColormap` only if it
  has fewer than :rcraw:`cmap.listedthresh` levels (:commit:`ade787f9`). This is
  critical in case users import cmaps from other projects.
* Permit constructing property cycles with `~proplot.constructor.Cycle` by passing
  ``color`` as keyword argument (:commit:`86a50eb2`). This is matplotlib-like workflow.
* Permit disabling property cycling with e.g. ``cycle=False``, ``cycle='none'``,
  or ``cycle=()``, and re-enabling the default with ``cycle=True`` (:commit:`86a50eb2`).
* Override `~matplotlib.axes.Axes.set_prop_cycle` to pass the input arguments
  through `~proplot.constructor.Cycle` (:commit:`86a50eb2`). Features are a superset
  and this also lets me cache the cycler for comparison with on-the-fly inputs.
* Add shorthands :rcraw:`grid.width`, :rcraw:`grid.style`, :rcraw:`gridminor.width`,
  and :rcraw:`gridminor.style` for the respective ``linewidth`` and ``linestyle``
  settings (:commit:`a50d5264`)
* Permit "registering stuff" by passing files or objects to
  `~proplot.config.register_cmaps`, `~proplot.config.register_cycles`,
  `~proplot.config.register_colors`, and `~proplot.config.register_fonts`
  rather than forcing users to use the ``.proplot`` folder (:commit:`ad999e95`).
* Support case insensitivity when calling matplotlib's ``unregister_cmap``
  by improving `~proplot.colors.ColormapDatabase` so it derives from a
  `~collections.abc.MutableMapping` rather than `dict` (:commit:`ade787f9`).
* Make `~proplot.config.Configurator` a `~collections.abc.MutableMapping`
  subclass and change iteration behavior to loop over keys, not item
  pairs (:commit:`5626bc88`).
* Rename `~proplot.config.Configurator.get` to `~proplot.config.Configurator.find`
  (:commit:`e8559f3d`). Confusing since ``get`` didn't accept a "fallback" second
  positional argument. Now ``get`` is the "dictionary-like" inherited method.
* Add public `~proplot.config.Configurator.changed` property to display a dictionary
  of settings changed from proplot defaults (:commit:`41b5e400`).
* Add public `~proplot.config.Configurator.user_file` and
  `~proplot.config.Configurator.user_folder` static methods for displaying
  folder locations (:commit:`b11d744a`).
* Support XDG directories for proplot config files on Linux (:issue:`204`,
  :commit:`5e6367dc`). Also accept the file ``~/.proplotrc`` and the folder
  ``~/.proplot`` on all systems and raise a warning if duplicate valid files
  or folders are found.
* Make `~proplot.config.rc_proplot` and `~proplot.config.rc_matplotlib` containers
  of proplot/matplotlib settings part of the public API (:commit:`a50d5264`).
* Allow conversion of numeric inputs with `~proplot.utils.units` using e.g.
  ``pplt.units(num, 'in', 'cm')`` (:commit:`88f3dc88`).
* Add more intuitive :rcraw:`grid.labelpad` and :rcraw:`tick.labelpad`
  as aliases for :rcraw:`grid.pad` and :rcraw:`tick.pad` (:commit:`a50d5264`).
* Add `~proplot.axes.PlotAxes.line` and `~proplot.axes.PlotAxes.linex` command
  aliases for `~proplot.axes.PlotAxes.plot` and `~proplot.axes.PlotAxes.plotx`
  (:commit:`6382cf91`). This is more intuitive.
* Add `~proplot.axes.PlotAxes.stepx` and `~proplot.axes.PlotAxes.stemx` commands
  analogous to `~proplot.axes.PlotAxes.plotx`, and add `~proplot.axes.PlotAxes.histh`,
  `~proplot.axes.PlotAxes.boxploth` (shorthand `~proplot.axes.PlotAxes.boxh`),
  and `~proplot.axes.PlotAxes.violinploth` (shorthand `~proplot.axes.PlotAxes.violinh`)
  commands analogous to `~proplot.axes.PlotAxes.barh` (:commit:`6382cf91`).
* Exclude out-of-bounds data when determining automatic y (x) axis limits when x (y)
  limits have been explicitly set for `plot` and `scatter` plots (:commit:`6382cf91`).
  Controlled by the :rcraw:`axes.inbounds` property, analogous to :rcraw:`cmap.inbounds`
  used for cmap scaling. This feature leverages proplot's input standardization.
* Let 1D plotting commands iterate over columns of 2D *x* and *y* coordinate arrays
  instead of only 2D *y* coordinate arrays (:commit:`6382cf91`.)
* Only modify `heatmap` major and minor tick locations if the
  default tickers are active (:pr:`6382cf91`).
* Support more artist synonyms throughout plotting overrides, e.g. ``ec``
  for `edgecolor`, `lw` for `linewidth`, `fc` and `fillcolor` for
  `facecolor` (:commit:`6382cf91`). This expands matplotlib synonyms.
* Support passing positional fifth-argument colors to `~proplot.axes.PlotAxes.barbs`
  and `~proplot.axes.PlotAxes.quiver`, just like `~proplot.axes.PlotAxes.scatter`
  (:commit:`6382cf91`). This was previously not possible.
* Use default luminance of ``90`` rather than ``100`` for auto-colormaps generated
  for barb, scatter, and streamline plots (:commit:`6382cf91`).
* Support automatic labels for ``tricontour`` and ``tripcolor`` plots alongside
  the more common ``contour`` and ``pcolor``. (:commit:`6382cf91`).
* Add `rasterize` keyword to `colorbar` so that colorbar solids rasterization can
  be turned on (proplot turns off by default) (:commit:`6382cf91`).
* Add `edgefix` keyword to `colorbar` to control colorbar-solid edges and
  use shared ``_fix_edges`` function (:commit:`6382cf91`).
* Add `location` keyword as alternative to `loc` for legend and
  colorbar funcs (:commit:`5cb839fd`).
* Add `alphabetize` keyword to `legend` to optionally alphabetize handles by
  their labels (:commit:`6382cf91`).
* Apply auto-detected xarray and pandas legend/colorbar titles even if the
  legend/colorbar are not drawn on-the-fly (:issue:`6382cf91`).
* Add :rcraw:`colorbar.facecolor` and :rcraw:`colorbar.edgecolor` properties
  analogous to legend properties for controlling frame (:commit:`6382cf91`).
* Treat singleton lists and tuple `legend` input same as scalar
  handle input, i.e. never triggers "centered row" specification (:commit:`6382cf91`).
* Support auto-detection of tuple-grouped `legend` handle labels when labels
  not passed explicitly (:commit:`6382cf91`).
* Automatically pull out grouped tuples of artists passed to `legend` if they have
  differing labels (:commit:`6382cf91`). This is useful for passing error shading to `legend`.
* Silently ignore non-artist and non-container `legend` input -- e.g., ignore the bins
  and values returned by `hist` (:commit:`6382cf91`).
* Allow list-of-list "centered row" `legend` specification with e.g.
  ``[h, [h1, h2, h3]]`` (i.e., mixed list and non-list input) (:commit:`6382cf91`).
* Permit partial specification of `legend` labels, e.g. ``[h1, h2]`` paired
  with ``['label', None]`` overrides the artist label for ``h1`` but uses
  the artist label for ``h2`` (:commit:`6382cf91`).
* Interpret all native matplotlib `legend` spacing arguments (e.g., `borderpad`
  and `columnspacing`) with `~proplot.utils.units` (:commit:`6382cf91`).
* Make `proplot.axes.Axes.colorbar` capture matplotlib-native `format`
  keyword as alias for `formatter` and `ticklabels` (:issue:`262`).
* Support list-of-string parametric coordinates and format on-the-fly colorbar
  ticks with those string labels (:commit:`02fbda45`). This may be a common
  use case for parametric plots.
* Add `ignore` keyword to omit specific ``show_cmaps``, ``show_cycles``, and
  ``show_colors`` categories from the tables (:issue:`c45d5fa1`).
* Allow case-insensitive specification of ``show_cmaps``, ``show_cycles``, and
  ``show_colors`` categories and never ignore input colormaps even if they
  match an ignored name like ``'jet'`` (:issue:`c45d5fa1`).
* Support restricting cartopy bounds in cartopy 0.19 by leveraging the
  `ylim` `~cartopy.mpl.gridliner.Gridliner` property (:commit:`e190b66c`).
* Add `xlabelpad`, `ylabelpad`, `xticklabelpad`, `yticklabelpad` keywords
  to `~proplot.axes.CartesianAxes.format` and read and apply changed
  :rcraw:`axes.labelpad` (:commit:`e7d86b8f`).
* Add support for "minor" radial and azimuthal gridlines in
  `proplot.axes.PolarAxes.format` (:commit:`59c85f0e`).
* Add `thetagrid`, `rgrid`, `thetagridminor`, and `rgridminor` keys to
  `proplot.axes.PolarAxes.format` to toggle gridlines, and read and apply changed
  toggles from rc settings -- consistent with Cartesian axes (:commit:`59c85f0e`).
* Add `title_kw`, `suptitle_kw`, `leftlabels_kw`, `rightlabels_kw`, `toplabels_kw`,
  and `bottomlabels_kw` to `proplot.axes.Axes.format` for arbitrarily modifying
  label text objects -- consistent with `xlabel_kw` and `ylabel_kw` used
  for `proplot.axes.CartesianAxes.format` (:commit:`6382cf91`).

Bug fixes
---------

* Fix issue with unpacking iterables inside return statements in python < 3.8
  (:pr:`268`) by `Eli Knaap`_.
* Fix issue where auto layout algorithm recurses in popup backends (:commit:`51967ce3`).
* Fix issue where auto layout algorithm blows up in mpl 3.4+ (:commit:`51967ce3`).
* Fix issue where tight layout is effectively deactivated in mpl >= 3.4 due to
  ``set_position`` automatically calling ``set_in_layout(False)`` (:commit:`20502345`).
* Fix issue where thin pyplot-function wrappers e.g. ``isinteractive``
  do not return results (:commit:`e62e3655`).
* Fix issue where `proplot.config.Configurator.save` preserves the ``'#'``
  in HEX strings, resulting in values that cannot be read back in with
  `proplot.config.Configurator.load` (:commit:`41b5e400`).
* Fix issue where deprecated `aspect` `~proplot.ui.subplots` argument
  is ignored (:commit:`70a8b87d`).
* Fix issue where explicit user-input ``width`` is ignored when creating
  colorbars or panels and gridspec slot already exists (:commit:`51967ce3`).
* Fix bug where the default space selection failed to use the
  figure-wide share setting (:commit:`51967ce3`).
* Fix bug where the reference subplot aspect ratio not preserved in
  presence of complex geometry with panels (:commit:`51967ce3`).
* Fix issue where a-b-c labels are removed in presence of ``'top'`` panels
  with ``titleabove=True`` (:commit:`7873d5e0`).
* Fix issue where 'aligned' labels fail in recent matplotlib versions
  due to private matplotlib API change (:commit:`51967ce3`).
* Fix issue where ``cmap.reverse()`` returns strange monochrome colormaps
  when channel values are specified by functions (e.g., ``cubehelix``) due
  to loop scope overwriting a non-local lambda function variable (:commit:`ade787f9`).
* Fix issue where ``_restrict_inbounds`` fails for reversed/descending axis
  limits (:commit:`6382cf91`).
* Fix issues where cartopy minor gridlines are toggled on when map bounds are changed
  and basemap map boundary props cannot be modified (:commit:`c1f1a7de`).
* Turn off ``_restrict_inbounds`` for geographic projections to prevent issue where
  lon/lat coordinates are compared to map coordinates (:commit:`6382cf91`). In-bounds
  colormap scaling for geographic projections may be added in a future version.
* Fix issue where error indications do not ignore masked values
  in masked numpy arrays (:commit:`6382cf91`).
* Fix issue where error shading objects are grouped into lists rather than tuples
  and are not combined into single handle when passed to ``legend`` (:issue:`260`).
* Fix issue where `~proplot.axes.Axes.parametric` ignores `interp` when
  selecting `DiscreteNorm` colormap levels (:commit:`152a3a81`).
* Fix issue where tight layout padding is not respected for panels created from
  twin axes by ensuring panel parent is always the main axes (:commit:`e7d86b8f`).
* Fix obscure bug where axis labels in presence of mixed panels and
  non-panels are improperly shared (:commit:`06666296`).
* Stop overwriting user-input `spineloc` when combined with user-input
  spine `bounds` (:commit:`e7d86b8f`).
* Include *children* of ``key`` when triggering complex synced settings
  (e.g., now we trigger application of :rcraw:`tick.widthratio` when either
  :rcraw:`tick.width` or :rcraw:`meta.width` are changed) (:commit:`5626bc88`).

Internals
---------

* Convert all plotting wrappers to dedicated overrides of individual functions
  in `~proplot.axes.PlotAxes` class (:commit:`6382cf91`). This massively simplifies
  the internals and makes learning and adopting proplot much easier for users.
* Implement interpretation of physical units and "panel" accounting directly
  on the `~proplot.gridspec.GridSpec` instead of cumbersome hidden methods
  in `~proplot.figure.Figure` (:commit:`20502345`).
* Validate all setting assignments to `~proplot.config.Configurator` using a new
  `~proplot.config.rc_proplot` dictionary, analogous to ``rcParams``
  (:pr:`109`, :commit:`5626bc88`). This helps prevent bugs.
* Move ``text``, ``legend``, and ``colorbar`` overrides to base `~proplot.axes.Axes`
  class separate from `~proplot.axes.PlotAxes` (:commit:`6382cf91`).
* Automatically redirect all internal plotting calls to native matplotlib methods
  (:commit:`6382cf91`). This significantly improves stability.
* Move ``register_colors`` internals from ``config.py`` to ``colors.py``
  by breaking up into smaller functions (:commit:`ad999e95`).
* Move ``_version`` to a separate ``dependencies.py`` file and
  allow more versatile comparison operations (:commit:`8806631d`).
* Efficiently impose `~proplot.axes.GeoAxes` defaults ``latlon=True`` and
  ``transform=PlateCarree()`` in 90% fewer lines by looping over funcs.

Documentation
-------------

* Move all plotting wrapper documentation to dedicated methods and remove
  references to wrappers in User Guide and Getting Started.
* Embed `proplot.figure.Figure` documentation inside `proplot.ui.subplots`
  instead of just referencing it.
* Embed `proplot.axes.Axes.format` documentation inside ``format``
  documentation for subclasses instead of just referencing it.
* Document the relative font size scalings with a table in
  `~proplot.axes.Axes.text` (:commit:`6382cf91`).

ProPlot v0.7.0 (2021-07-11)
===========================

Deprecated
----------

* Remove v0.6.0 renamed classes (e.g. `ProjAxes`) from top-level namespace
  (:commit:`442e6aa6`). These were kept available just for documentation. The renamed
  functions `shade`, `saturate`, and `inline_backend_fmt` remain available.
* Change default :rcraw:`savefig.transparent` back to ``False`` (:pr:`252`). Dubious
  justification for ``True`` in the first place, and makes default PNG proplot figures
  unreadable wherever "dark mode" is enabled.
* Rename SciVisColor colormaps from ``Blue1``, ``Blue2``, etc. to plurals ``Blues1``,
  ``Blues2``, etc. to avoid name conflict with open-color colors (:commit:`8be0473f`).
  Requesting the old names (case-sensitive) redirects to the new names
  (:commit:`3f0794d0`). This permits making monochromatic open-color maps with e.g.
  ``plot.Colormap('blue9')`` and feels more consistent with ColorBrewer convention of
  using plurals like ``Blues``, ``Reds``, etc.
* Shuffle various SciVisColor colormap names to make them consistent/succinct. Make
  ``Browns1`` the most colorful/vibrant one, just like ``Greens1`` and ``Blues1``;
  split up the ``RedPurple`` maps into ``Reds`` and ``Purples``; and add
  the ``Yellows`` category from the ``Oranges`` maps (:commit:`8be0473f`). Requesting
  the old names (case-sensitive) redirects to the new names (:commit:`3f0794d0`).
* Add :rcraw:`image.discrete` options and `discrete` keyword for toggling
  `~proplot.colors.DiscreteNorm` application, and disable by default for `imshow`,
  `matshow`, `spy`, `hexbin`, and `hist2d` plots (:issue:`233`, :commit:`5a7e05e4`).
  Also make `hexbin` and `hist2d` behavior with ``discrete=True`` more sane by using
  maximum possible counts for autoscaling, and change `~proplot.colors.DiscreteNorm`
  argument `extend` to more intuitive name `unique`.
* Rename :rcraw:`subplots.pad` and :rcraw:`subplots.axpad` to more intuitive
  :rcraw:`subplots.outerpad` and :rcraw:`subplots.innerpad` (:commit:`3c7a33a8`).
  Also rename `~proplot.figure.Figure` keywords.
* Rename `width` and `height` `~proplot.subplots.subplots` keyword args to `figwidth`
  and `figheight` to avoid confusion with `refwidth`/`refheight` (:commit:`12d01996`).
  Will accept old keyword args without warning since they are used heavily.
* Rename `aspect`, `axwidth`, and `axheight` keyword args to more intuitive
  `refaspect`, `refwidth`, and `refheight` (:commit:`12d01996`). Will accept old
  keyword args without warning since they are used heavily.
* Rename `abovetop` keyword for moving title/abc labels above top panels, colorbars,
  and legends to :rcraw:`title.above` (:commit:`9ceacb7b`). Example usage:
  ``ax.format(title='Title', titleabove=True)``.
* Rename the `proplot.colors.PerceptuallyUniformColormap.from_color` keywords `shade`,
  `fade` to `luminance`, `saturation` keyword (:commit:`3d8e7dd0`). These can also
  be passed to `~proplot.contructor.Colormap` when it is called with positional arguments.
* Rename seldom-used `Figure` argument `fallback_to_cm` to more understandable
  `mathtext_fallback` (:pr:`251`).
* Reduce default :rcraw:`savefig.dpi` to 1000 (:commit:`bfda9c98`). Nature recommends
  1000, Science recommends "more than 300", PNAS recommends 1000--1200. So 1000 is fine.
* Increase default :rcraw:`colorbar.insetpad` to avoid recurring issue where ticklabels
  run close to the background patch (:commit:`f5435976`)
* Use proplot TeX Gyre fonts with `~proplot.config.use_style` styles unless specified
  otherwise (:commit:`6d7444fe`). Styles otherwise build on matplotlib defaults.
* When using ``medians=True`` or ``means=True`` with `indicate_error` plot simple
  error bars by default instead of bars and "boxes" (:commit:`4e30f415`). Only plot
  "boxes" with central "markers" by default for violin plots (:commit:`13b45ccd`).
* `legend_extras` no longer returns the background patch generated for centered-row
  legends (:pr:`254`). This is consistent with `colorbar_extras` not returning
  background patches generated for inset colorbars. Until proplot adds new subclasses,
  it makes more sense if these functions only return `~matplotlib.legend.Legend` and
  `~matplotlib.colorbar.Colorbar` instances.

Features
--------

* Add the remaining commonly-used backend-related `pyplot` functions `ion`, `ioff`,
  `isinteractive`, and `switch_backend` to the top-level `proplot` namespace
  (:commit:`cd440155`). This avoids forcing users to import pyplot inside a proplot
  session (the remaining pyplot functions are related to the "non-object-oriented"
  workflow, which proplot explicitly discourages).
* Add support for local ``proplotrc`` files in addition to "hidden"
  ``.proplotrc`` files with leading dot (:commit:`8a989aca`).
* Add minimal support for "3D" `~matplotlib.mpl_toolkits.mplot3d.Axes3D` axes
  (:issue:`249`). Example usage: ``fig.subplots(proj='3d')``.
* Add `wequal`, `hequal`, and `equal` options to still use automatic spacing but force
  the tight layout algorithm to make spacings equal (:pr:`215`, :issue:`64`)
  by `Zachary Moon`_.
* Determine colormap levels using only in-bounds data if the *x* or *y* axis limits
  were explicitly set (:issue:`209`). Add `inbounds` `~proplot.axes.apply_cmap`
  keyword and :rcraw:`image.inbounds` setting to control this.
* Allow calling `proplot.colors.PerceptuallyUniformColormap.from_hsl` by passing
  `hue`, `saturation`, or `luminance` to `~proplot.constructor.Colormap` without
  any positional arguments (:commit:`3d8e7dd0`).
* Allow passing `alpha`, `luminance`, `saturation` to `~proplot.constructor.Colormap`
  as lists to be applied to each component cmap (:commit:`3d8e7dd0`).
* Add convenient shorthands for channel references throughout colormap functions --
  e.g. `h` for hue, `l` for `luminance`, etc. (:commit:`3d8e7dd0`).
* Add the ``'Flare'`` and ``'Crest'`` seaborn colormaps (:commit:`14bc16c9`). These
  are seaborn's color cycle-friendly alternatives to existing maps.
* Add the `~proplot.utils.shift_hue` function analogous to `scale_saturation`
  and `scale_luminance` (:commit:`67488bb1`).
* Add the `~proplot.utils.to_hex` function and make all color-manipulation funcs return
  HEX strings by default (:commit:`67488bb1`). Otherwise `scatter` throws warnings.
* Use ``90`` as the default `luminance` when creating monochromatic colormaps with
  `to_listed` set to ``True`` (as when `~proplot.constructor.Cycle` calls
  `~proplot.constructor.Colormap`; :commit:`3d8e7dd0`).
* Add `~proplot.axes.Axes.plotx` and `~proplot.axes.Axes.scatterx` commands that
  interpret plotting args as ``(y, x)`` rather than ``(x, y)``, analogous to
  `~proplot.axes.Axes.areax` (:pr:`258`).
* Add support for `~proplot.axes.indicate_error` *horizontal* error bars and shading
  for *horizontal* plotting commands `barh`, `plotx`, and `scatterx` (:pr:`258`).
* Add support for ``ax.plot_command('x_key', 'y_key', data=dataset)`` for
  virtually all plotting commands using `standardize_1d` and `standardize_2d`
  (:pr:`258`). This was an existing `~matplotlib.axes.Axes.plot` feature.
* Add support for the plotting style ``ax.plot(x1, y1, fmt1, x2, y2, fmt2, ...)``
  as allowed by matplotlib (:pr:`258`).
* Add `absolute_width` keyword to `~proplot.plot.bar_extras` to make `width`
  argument absolute (:pr:`258`). Remains ``False`` by default.
* Use "sticky" edges in x-direction for lines drawn with `plot()` and in y-direction
  for lines drawn with `plotx()` (:pr:`258`). This eliminates padding along the
  "dependent" axis when limits are not specified, similar to histograms and
  barplots and matching a feature we previously added to `fill_between` (:pr:`166`).
* Add support for "stacked" plots to `~matplotlib.axes.Axes.vlines` and
  `~matplotlib.axes.Axes.hlines` (:pr:`258`).
* Add `stack` as alternative to `stacked` for bar and area plots (:commit:`4e30f415`).
  Imperative keywords are better.
* Allow passing e.g. ``barstds=3`` or ``barpctiles=90`` to request error bars
  denoting +/-3 standard deviations and 5-95 percentile range (:commit:`4e30f415`).
* Add singular `indicate_error` keywords `barstd`, `barpctile`, etc. as
  alternatives to `barstds`, `barpctiles`, etc. (:commit:`81151a58`).
  Also prefer them in the documentation.
* Permit different colors for `~matplotlib.axes.Axes.boxplot` and
  `~matplotlib.axes.Axes.violinplot` using color lists (:issue:`217`, :pr:`218`)
  by `Mickaël Lalande`_. Also allow passing other args as lists (:commit:`4e30f415`).
* Allow passing ``means=True`` to `boxplot` to toggle mean line
  (:commit:`4e30f415`).
* Allow setting the mean and median boxplot linestyle with
  ``(mean|median)(ls|linestyle)`` keywords (:commit:`4e30f415`).
* Automatically set ``fill=True`` when passing a fill color or color(s)
  to `boxplot_wrapper` (:commit:`4e30f415`).
* Allow updating `vlines` and `hlines` styling with singular `color` and `linestyle`
  and all of their aliases (:pr:`258`).
* Allow updating axes fonts that use scalings like ``'small'`` and ``'large'``
  by passing ``fontsize=N`` to `format` (:issue:`212`).
* Add `titlebbox` and `abcbbox` as alternatives to `titleborder` and `abcborder` for
  "inner" titles and a-b-c labels (:pr:`240`) by `Pratiman Patel`_. Borders are still
  used by default.
* Allow putting `title` and `abc` in the same location -- the title and label
  are simply offset away from ech other (:issue:`402214f9`). Padding between
  them is controlled by the new param :rcraw:`abc.titlepad`.
* Add new :rcraw:`suptitle.pad`, :rcraw:`leftlabel.pad`, :rcraw:`toplabel.pad`,
  :rcraw:`bottomlabel.pad`, :rcraw:`rightlabel.pad` settings to control padding
  used when aligning super labels (:commit:`402214f9`). These can also be passed
  to `~proplot.axes.Axes.format` and applied locally. The new defaults increase
  super title padding by a bit.
* More robust interpretation of :rcraw:`abc.style` -- now match case with first
  ``'a'`` or ``'A'`` in string, and only replace that one (:issue:`201`).
* Interpret fontsize-relative legend rc params like ``legend.borderpad``
  with ``'em'`` as default units rather than ``'pt'`` (:commit:`6d98fd44`).
* Add :rcraw:`basemap` setting for changing the default backend (:commit:`c9ca0bdd`). If
  users have a cartopy vs. basemap preference, they probably want to use it globally.
* Add :rcraw:`cartopy.circular` setting for optionally disabling the "circular bounds
  on polar projections" feature (:commit:`c9ca0bdd`).
* Support the standard aliases ``'ls'``, ``'linestyle'``, ``'linestyles'``, etc.
  in `~proplot.constructor.Cycle` calls (:commit:`3d8e7dd0`).
* Add `queue` keyword to `colorbar` and `legend` to support workflow where users
  successively add handles to location (:pr:`254`).
* Add `nozero` keyword arg to `apply_cmap` to remove the zero contour
  from automatically generated levels (:commit:`10e0f13b`).
  Example usage: ``ax.contour(x, y, z, nozero=True)``.
* Add `positive` and `negative` keyword args to `apply_cmap` for requesting
  automatically-generated all-positive or all-negative levels (:commit:`335d58f4`).
  Example usage: ``ax.contourf(x, y, z, positive=True)``.
* Add `rotation` keyword to `colorbar_wrapper` for rotating colorbar tick
  labels, like `xrotation` and `yrotation` (:commit:`2d835f20`).
* Add `tickdir` and `tickdirection` keywords to `colorbar_wrapper` for
  controlling tick style, like `xtickdir` and `ytickdir` (:commit:`f377f090`).
* Allow specifying labels for auto-generated legends using a ``'labels'`` key
  in a `legend_kw` keyword argument (:commit:`a11d1813`).
* Replace legends drawn in the same location by default rather than drawing two
  legends on top of each other (:pr:`254`).
* Use `Artist` labels for the default list-of-artist colorbar tick labels if `values`
  was not passed -- and if labels are non-numeric, rotate them 90 degrees for horizontal
  colorbars by default (:commit:`ed8e1314`). Makes the choice between "traditional"
  legends and "colorbar-style" legends more seamless.
* Use same default-level generation algorithm for contour plots without colormaps as for
  all other colormap plots (:commit:`10e0f13b`). Makes automatically-generated
  solid-color contours and colormap-style contours identical.
* Add suffix ``'_copy'`` to colormaps converted with `to_listed` and
  `to_linear_segmented` to avoid accidental overwriting (:commit:`91998e93`).
* If available, use :rcraw:`pcolormesh.snap` to repair overlap in transparent colorbar
  solids rather than manual-blending workaround (:commit:`c9f59e49`).
* Add `xmin`, `xmax`, `ymin`, and `ymax` keyword args to
  `~proplot.axes.CartesianAxes.format` as alternatives to `xlim` and `ylim`
  (:commit:`ae0719b7`). Example usage: ``ax.format(xmin=0)`` as opposed to
  ``ax.format(xlim=(0, None))``.
* Allow passing full "side" names to `lonlabels` and `latlabels` rather than
  abbreviations, e.g. ``'left'`` instead of ``'l'`` (:commit:`a5060f67`). This is
  more consistent with rest of package.
* Set default transform to ``ccrs.PlateCarree`` when calling `matplotlib.axes.Axes.fill`
  on `CartopyAxes` (:issue:`193`). This is more consistent with rest of package.

Bug fixes
---------

* Fix 3 fatal issues preventing proplot import and basic usage in matplotlib >= 3.4
  (:pr:`251`).
* Fix deprecation warnings associated with matplotlib 3.4 refactoring of
  subplot classes (:pr:`251`).
* Fix deprecated reference to :rc:`fallback_to_cm` in matplotlib >= 3.3
  (:pr:`251`).
* Fix `~matplotlib.ticker.IndexFormatter` deprecation warning in matplotlib >= 3.3 by
  replacing with proplot-local copy (:pr:`251`).
* Fix deprecation warning in matplotlib >= 3.3 -- add `extend` as mappable attribute
  rather than passing it to `colorbar()` (:commit:`a23e7043`).
* Fix issue where figures with fixed-aspect axes don't scale properly
  in matplotlib >= 3.3 (:issue:`210`, :issue:`235`).
* Fix issue where "twin" ("alternate") axes content always hidden beneath "parent"
  content due to adding as children (:issue:`223`).
* Fix issue where default layout in complex subplot grids with non-adjacent
  edges is incorrect (:issue:`221`).
* Fix issue where `apply_cycle` fails to merge mean-uncertainty legend handles
  due to presence of placeholder labels (:commit:`4e30f415`).
* Fix issue where `standardize_1d` inappropriately infers legend entries from
  y-coordinate metadata rather than column metadata (:commit:`4e30f415`).
* Fix issue where `barb` and `quiver` cannot accept 1D data arrays (:issue:`255`).
* Fix issue where cannot set ``rc.style = 'default'`` (:pr:`240`) by `Pratiman Patel`_.
* Fix issue where `get_legend` returns None even with legends present (:issue:`224`).
* Fix issue where new child axes reset row/col label settings (:commit:`f32d9703`).
* Fix issue where `~xarray.DataArray` string coordinates are not extracted from
  container before applying as tick labels (:issue:`214`).
* Fix issue where cannot set `extend` other than ``'neither'`` for
  `~matplotlib.axes.Axes.scatter` colorbars (:issue:`206`).
* Fix issue where `~matplotlib.axes.Axes.hexbin` ignores `vmin` and `vmax`
  keywords (:issue:`250`).
* Fix issue where parametric plot *x* axis is reversed (:commit:`3bde6c47`).
* Fix issue where e.g. `ax.area(x, 0, y2, negpos=True` has positive colors
  below x-axis and negative above x-axis (:pr:`258`).
* Fix issue where "negpos" plots ignore `edgecolor` because they pass
  `color` rather than `facecolor` to plotting commands.
* Fix issue where cannot have datetime labels on `area` plots (:issue:`255`).
* Fix issue where default orientation of `barh` vertical axis is reversed
  (:commit:`258`).
* Fix issue where `hist` with `xarray.DataArray` or `pandas.Dataframe` input causes
  erroneous axis labels; use labels for legend instead (:issue:`195`).
* Fix issue where axis is accidentally inverted for histogram plots (:issue:`191`).
* Fix issue where `[xy]minorlocator=1` is not allowed (:issue:`219`).
* Fix issue where inner titles ignore axes-local `titlepad` (:commit:`14f3d0e3`).
* Fix issue where we again fail to sufficiently pad title above tick marks
  with tick marks on top x-axis (:commit:`402214f9`).
* Fix issue where non-Cartesian `heatmap` errors rather than warns (:issue:`238`).
* Fix issue where ``labels=True`` with no contours causes error (:issue:`238`).
* Fix issue where `~proplot.colors.Cycle` fails to register new names and fails to
  display in `~proplot.demos.show_cycles` (:commit:`94ffc1dc`, :commit:`4a7a3c79`).
* Fix issue where proplot ignores `set_under` and `set_over` values when translating
  matplotlib colormap classes to proplot subclasses (:issue:`190`).
* Fix issue where `~proplot.colors.DiscreteNorm` does not account for `set_under` and
  `set_over` colors distinct from adjacent in-bounds colors (:issue:`190`).
* Fix issue where proplot fails to detect legend entries for "outer"
  legends (:issue:`189`).
* Fix issue where list-of-list-style `legend()` handle and label input fails completely
  (:commit:`a298f81f`). This input style is used to specify "centered" legend rows.
* Fix error message when no legend handles are found (:commit:`2c6bf3e2`).
* Fix issue where multiple-artist legend entries (e.g., for lines indicating means and
  shading indicating uncertainty) are accidentally truncated (:commit:`a11d1813`).
* Fix issue where numeric zero cannot be applied as legend label (:commit:`02417c8c`).
* Fix issue where simple `pandas.DataFrame.plot` calls with ``legend=True`` fail
  (:pr:`254`, :issue:`198`).
* Fix unnecessary restriction where users can only draw <2 "alt" axes and clean
  up the `alt[xy]` and `dual[xy]` internals (:issue:`226`).
* Fix matplotlib bug where `altx` and `alty` reset the minor locator of the shared
  axis to ``AutoMinorLocator`` even if the axis scale is ``'log'`` (:commit:`2f64361d`).
* Fix issue where axis coordinates are incorrect when `violinplot` or `boxplot`
  receive non-DataFrame input (:commit:`b5c3ec4c`).
* Fix issue where `indicate_error` cannot accept 1D error bounds (:commit:`ef2d72cd`).
* Fix issue where `show_cmaps` cannot display reversed colormaps (:commit:`2dd51177`).
* Fix issue where ``'grays_r'`` translated to ``'greys'`` (:commit:`074c6aef`).
* First reverse, *then* shift ``cmap_r_s`` colormaps (:commit:`e5156294`).
* Fix obscure `~proplot.axes.Axes.parametric` bug where `numpy.stack` tries to make
  nested ragged arrays from parametric coords (:commit:`b16d56a8`).
* Fix issue where where `SubplotSpec.get_active_rows_columns` returned incorrect
  number of "active" rows and columns (:commit:`5cf20b84`).
* For rc lookup with ``context=True``, use most restrictive search mode rather than least.
  Otherwise `ax.format()` calls inside context blocks can be overwritten with the
  default rc values in subsequent `ax.format()` calls (:commit:`8005fcc1`).

Internals
---------

* Refactor massive `standardize_(1d|2d)` and `(cmap|cycle)_changer` wrappers to break
  things into manageable chunks (:pr:`258`, :commit:`6af22567`, :commit:`d3352720`).
* Refactor `colorbar` and `legend` methods and their massive wrappers to clean
  things up and expand the "queueing" feature beyond wrappers (:pr:`254`).
* Add prefix ``'proplot_'`` to registered axes "projections" (:commit:`be7ef21e`). More
  clear and guards against conflicts with external packages and other mpl versions.
* Add system for processing flexible keyword arguments across different commands
  to ``internals/__init__.py``. Analogous to mpl ``_alias`` processing.

Documentation
-------------

* Finally use ``pplt`` as the recommended abbreviation: ``import proplot as pplt``.
* Major clean up of "Why ProPlot?" page and user guide pages.
* Fix incomplete ``cmap.from_file`` docstrings (:commit:`54f1bc7c`).
* Rename "Changelog" to "What's New?" and list all contributors in "About the Authors".
* Rename public/documented funcs ending in `_wrapper` to ending in `_extras` to avoid
  implication they are the only funcs wrapping those commands (:commit:`d1e1e85b`).
* Rename public/documented func `make_mapping_array` to private function,
  following lead of matplotlib's `makeMappingArray` (:commit:`66ae574b`).
* Rename public/documented funcs `cmap_changer` and `cycle_changer`
  to `apply_cmap` and `apply_cycle` (:commit:`86f7699a`).


ProPlot v0.6.4 (2020-06-13)
===========================

Features
--------

* Change ``autoformat`` from a `Figure` keyword argument into the
  :rcraw:`autoformat` rc setting (:commit:`3a7e5a7c`).
* Combine shading and lines when drawing on-the-fly legends with `indicate_error`
  shading using tuple of `fill_between`, `plot` handles, and have `shadelabel` and
  `fadelabel` instead create separate entries *only when passed* (:issue:`187`).

Bug fixes
---------

* Fix major issue where calling ``legend()`` without any handles
  triggers error rather than using default handles (:issue:`188`).
* Fix issue where on-the-fly colorbar labels were
  ignored (:commit:`a642eeed`).
* Stop overwriting existing axis labels when ``autoformat=True``
  and DataArrays or DataFrames passed to plotting command (:commit:`76c7c586`).
* Support single-level contours with colormap colors (:issue:`182`).
* Support changing line width, line style, and color properties
  for barb, quiver, streamplot, matshow, spy, and hist2d plots
  (:issue:`177`).
* Use :rcraw:`patch.linewidth` for default bar edge width, stop setting
  default histogram plot linewidth to zero, and set :rcraw:`patch.linewidth`
  to ``0.6`` to match proplot's default line width for lines, axes edges, and
  hatches (:issue:`186`).

ProPlot v0.6.3 (2020-06-02)
===========================

Bug fixes
---------

* Fix issue where proplot import fails if cartopy is not installed (:commit:`e29d49e8`).

ProPlot v0.6.2 (2020-06-02)
===========================

Deprecated
----------

* Remove `~proplot.figure.Figure` setters like `set_sharex`, replace with
  read-only properties (:commit:`7b455008`). These did not work and did not
  add critical functionality.

Features
--------

* Add `autoformat` as `~proplot.axes.standardize_1d` and
  `~proplot.axes.standardize_2d` keyword arg, so inheriting labels can
  be turned on/off for individual plots (:commit:`61258280`).
* Share *initial* limits/scales/tickers from parent subplots when making
  new panels (:commit:`cf0d5d4e`).
* Permit negative "cuts" with `~proplot.colors.LinearSegmentedColormap.cut`
  to expand the neutral zone of a diverging cmap (:commit:`94548d09`).
* Add valid `format` arguments to `altx` and `alty`, including ``[x|y]lim``
  (:commit:`734f5940`).
* Pass string `dual[x|y]` arguments like ``'inverse'`` through the
  `~proplot.constructor.Scale` constructor (:commit:`413e1781`).
* Add ``'dms'`` locator and formatter, for degree-minute-second labels
  without cardinal direction indicators (:commit:`1b180cd2`).
* Add `"tau" formatter <https://tauday.com/tau-manifesto>`__
  (:commit:`fc6a9752`).
* Restore default :rcraw:`title.pad` to matplotlib value, stop artificially bumping
  up :rcraw:`title.pad` for "inner" titles (:commit:`7de1c1f4`).
* Make custom formatters like ``SciFormatter`` *classes* rather than functions
  returning `~matplotlib.ticker.FuncFormatter` (:commit:`7591f474`).

Bug fixes
---------

* Various improvements to auto-figure sizing with Qt backend and when calling
  `print_figure` (:commit:`db4e48d5`, :commit:`82457347`, :commit:`744d7d37`).
* Suppress warning when ``matplotlibrc`` contains non-style param
  (:commit:`4a0c7f10`).
* Fix fatal `standardize_2d` error when ``autoformat=False`` (:issue:`181`)
* Fix issue where ``Colormap(..., alpha=alpha)`` made persistent changes
  to the original registered colormap (:commit:`cb24ea51`).
* Prevent matplotlib deprecation warning by removing `set_smart_bounds`
  dependency and improving axis scale transforms (:commit:`432576d8`).
* Fix panel sharing issue in presence of stacked or multiple panels
  (:commit:`28eaf0ca`).
* Fix geographic feature toggling, zorder bugs (:commit:`acf0d5d4`, :commit:`ea151b25`).
* Fix `~matplotlib.axes.Axes.hist` bug due to ``bar(..., width=width)`` now
  being *relative* to the *x* step size (:commit:`e32ed0bc`).
* Fix bug where `~matplotlib.figure.Figure.savefig` receives ``Path`` instead
  of string (:issue:`176`).

Documentation
-------------

* Various improvements to website and API docstrings.
* Document `proplot.figure.Figure.save` method (:commit:`da25266a`).
* Darker "dark mode" (:commit:`979c8188`).
* Prevent website from flashing light mode when changing pages (:commit:`75e4d6a1`).

ProPlot v0.6.1 (2020-05-20)
===========================

Bug fixes
---------

* Fix issue where cartopy version checking fails if cartopy is not installed
  (:commit:`86cd50b8`).
* Fix issue where "tight" layout of geographic plots was broken in pre-v0.18
  cartopy (:commit:`72cb93c6`).
* Fix issue where gridline coverage was incomplete in some zoomed-in
  projections (:commit:`458c6d7c`).
* Fix issue where basemap minor gridlines did not update when
  major gridlines were updated (:commit:`427326a7`).

ProPlot v0.6.0 (2020-05-20)
===========================

Deprecated
----------

* Remove the ``geoaxes`` and ``geogrid`` rc settings (:pr:`168`). Gridline
  settings are now controlled with ``grid``.
* Remove the ``lonstep`` and ``latstep`` settings -- we now use
  `~proplot.ticker.LongitudeLocator` and `~proplot.ticker.LatitudeLocator`
  to select "nice" gridline locations even when zoomed in (:pr:`168`)
* Change default rc settings closer to matplotlib, including margins and line
  width (:pr:`166`, :commit:`f801852b`). Many were changed for no good reason.
* Change default line style for geographic gridlines from ``':'`` to ``'-'``
  and match style from primary gridlines (:pr:`166`, :commit:`f801852b`).
* Rename `add_errorbars` to `~proplot.axes.indicate_error` and rename
  various keyword args (:pr:`166`, :commit:`d8c50a8d`).
* Remove ``'rgbcycle'`` setting (:pr:`166`, :commit:`6653b7f0`).
* Deprecate support for "parametric" plots inside `~matplotlib.axes.Axes.plot`,
  instead use `~proplot.axes.Axes.parametric` (:commit:`64210bce`).
* Change `~proplot.utils.units` ``units`` keyword argument to more natural
  ``dest`` (:commit:`62903b48`).
* Remove the public objects `normalizers`, `locators`, `formatters`,
  `cartopy_projs`, `basemap_kwargs`, `cmaps`, `colors`, and `fonts` (:pr:`149`).
* Drop support for ``.xrgb`` and ``.xrgba`` files (:commit:`4fa72b0c`).  Not
  sure if any online sources produce these kinds of files.
* Drop support for ``.rgba`` files, but optionally read 4th opacity column
  from ``.rgb`` and ``.txt`` files (:commit:`4fa72b0c`).
* Stop reversing the ``'Spectral'`` colormap when ProPlot is imported
  (:pr:`149`, :commit:`ce4ef6a0`).
* Remove ``'Blue0'`` SciVisColor colormap (:pr:`149`, :commit:`7cb4ce0f`). It was odd
  man out in the table, and not even really perceptually uniform.
* Remove custom ProPlot cycles -- these should be thought out much more
  carefully (:commit:`43f65d17`).
* Remove "crayola" colors and clean up the `~proplot.setup.register_colors` algorithm
  (:pr:`149`, :commit:`8922d6de`). Crayola color names less intuitive than XKCD.
* Use ``'cmap_s'`` instead of ``'cmap_shifted'`` to quickly get a 180
  degree-shifted colormap, similar to ``'cmap_r'`` (:pr:`149`, :commit:`da4ccb08`).
* Rename ``GrayCycle`` colormap to ``MonoCycle`` to more accurately reflect
  colormap design origins (:pr:`149`, :commit:`d67e45bf`).
* Rename `~proplot.config.rc_configurator` and `~proplot.ui.subplot_grid` to
  `~proplot.config.RcConfigurator` and `~proplot.ui.SubplotsContainer`
  to match capitalized class naming convention (:pr:`149`).
* Rename `~proplot.colors.MidpointNorm` to more intuitive
  `~proplot.colors.DivergingNorm`, and make "fair" color scaling the default
  behavior (:commit:`2f549c9`).
* Rename `XYAxes` to `~proplot.axes.CartesianAxes`, `~proplot.axes.GeoAxes`
  to `~proplot.axes.CartopyAxes`, and `~proplot.axes.ProjAxes` to
  `~proplot.axes.GeoAxes` (:pr:`149`, :commit:`4a6a0e34`).
* Rename `BinNorm` to `~proplot.styletools.DiscreteNorm`
  and fix issues with diverging norm color scaling (:pr:`149`, :commit:`98a976f1`).
* Rename `ColorDict` to `~proplot.colors.ColorDatabase`, `CmapDict`
  to `~proplot.colors.ColormapDatabase` (:pr:`149`, :commit:`9d7fd3e0`).
* Rename `~proplot.styletools.LinearSegmentedColormap.concatenate` to
  `~proplot.styletools.LinearSegmentedColormap.append`,
  `~proplot.styletools.LinearSegmentedColormap.updated` to
  `~proplot.styletools.LinearSegmentedColormap.copy`,
  `~proplot.styletools.LinearSegmentedColormap.truncated` to
  `~proplot.styletools.LinearSegmentedColormap.truncate`, and
  `~proplot.styletools.LinearSegmentedColormap.punched` to
  `~proplot.styletools.LinearSegmentedColormap.cut` (:pr:`149`, :commit:`e1a08930`).
  The old method names remain with a deprecation warning.

Features
--------

* Add `~proplot.ticker.SigFigFormatter` (:pr:`149`, :commit:`da6105d2`)
  and `~proplot.ticker.SciFormatter` (:pr:`175`, :commit:`c43f7f91`)
  axis formatters.
* Make default `areax` and `areay` bounds "sticky", similar to
  histograms and barplots (:pr:`166`).
* Use `_LonAxis` and `_LatAxis` dummy axes with custom `LongitudeLocator`
  and `LatitudeLocator` to control geographic gridlines (:pr:`168`).
* Add ``'dmslat'`` and ``'dmslon'`` as formatters for cartopy projections,
  along with ``dms`` `format` keyword argument. This labels points with
  degrees/minutes/seconds when appropriate (:pr:`168`).
* Support "minor" geographic gridlines with the ``gridminor`` keyword
  arg and existing ``gridminor`` settings (:pr:`168`). Default locator
  used for minor gridlines is `~matplotlib.ticker.AutoMinorLocator`.
* Add `loninline`, `latinline`, and `rotatelabels` keywords for controlling
  cartopy gridliner behavior (:pr:`168`).
* Add `proplot.config.RcConfigurator.save` and
  `proplot.config.RcConfigurator.from_file` methods (:pr:`167`, :commit:`e6dd8314`).
* Increase default :rcraw:`savefig.dpi` to 1200, matching recommendations
  from academic journals (:pr:`167`, :commit:`c00e7314`). Also add detailed discussion
  to user guide.
* No longer distinguish between "quick" settings and proplot's "added"
  settings (:pr:`167`, :commit:`e6dd8314`). Quick settings, added settings, and
  matplotlib settings can all have "children" so the distinction no longer makes sense.
* Add opacity-preserving functions `~proplot.utils.to_rgba`
  and `~proplot.utils.to_xyza`, plus `~proplot.utils.set_alpha` for
  changing alpha channel of arbitrary color (:pr:`171`, :commit:`81c647da`).
* Add to `~proplot.colors.LinearSegmentedColormap.set_alpha` the ability to
  create an *opacity gradation*, rather than just an opacity for the entire
  colormap (:pr:`171`, :commit:`4583736`).
* Support passing colormap objects, not just names, to `~proplot.demos.show_cmaps`
  and `~proplot.demos.show_cycles` (:pr:`171`, :commit:`7f8ca59f`).
* Add options to `~proplot.axes.indicate_error` for adding *shading*
  to arbitrary plots (:pr:`166`, :commit:`d8c50a8d`). Also support automatic legend
  entries for shading and ensure `indicate_error` preserves metadata.
* Wrap ``pcolorfast`` just like ``pcolor`` and ``pcolormesh`` are
  wrapped (:pr:`166`, :commit:`50a262dd`).
* Add ``negpos`` feature to `~proplot.axes.bar_wrapper` and new :rcraw:`negcolor`
  and :rcraw:`poscolor` rc keyword arguments (:pr:`166`, :commit:`ab4d6746`).
* Support `~matplotlib.axes.Axes.vlines` and `~matplotlib.axes.Axes.hlines`
  flexible arguments and add ``negpos`` feature
  (:pr:`166`, :commit:`1c53e947`, :commit:`e42ee913`).
* Support `cartopy 0.18 <https://scitools.org.uk/cartopy/docs/latest/whats_new.html>`__
  locators, formatters, deprecations, and new labelling features (:pr:`158`).
* Support building a colormap and `DiscreteNorm` inside `~matplotlib.axes.Axes.scatter`,
  just like `contourf` and `pcolormesh` (:pr:`162`).
* Add :rcraw:`geogrid.labelpad` and :rcraw:`geogrid.rotatelabels` settings
  for cartopy gridline labels (:pr:`158`).
* Support more `~proplot.ticker.AutoFormatter` features on
  `~proplot.ticker.SimpleFormatter` (:pr:`152`, :commit:`6decf962`).
* Support drawing colorbars with descending levels (:pr:`149`, :commit:`10763146`)
* Add support for matplotlib stylesheets with `~proplot.config.use_style`
  function and ``style`` rc param (:pr:`149`, :commit:`edc6f3c9`).
* Add `categories` keyword arg to `~proplot.styletools.show_cmaps` and
  `~proplot.styletools.show_cycles` (:pr:`149`, :commit:`79be642d`).
* *Hide* bad colormaps like ``'jet'`` from the
  `~proplot.styletools.show_cmaps` table instead of deleting them outright,
  just like CSS4 colors (:pr:`149`, :commit:`ce4ef6a0`).
* Draw `~proplot.styletools.show_colors` table as single figure with category
  labels, similar to `~proplot.styletools.show_cmaps` (:pr:`149`, :commit:`c8ca2909`).
* Make ``'Grays'`` and ``'Greys'`` synonyms for the same ColorBrewer colormap
  (:pr:`149`, :commit:`da4ccb08`).
* Permit drawing "outer" axes and figure legends without explicitly passing
  handles (:pr:`149`, :commit:`a69b48eb`). Figure legends use the handles from all axes.
* Add `~proplot.styletools.LinearSegmentedColormap.to_listed` and
  `~proplot.styletools.PerceptuallyUniformColormap.to_linear_segmented`
  methods for handling conversions (:pr:`149`, :commit:`e1a08930`).
* Permit merging mixed colormap types `~proplot.styletools.LinearSegmentedColormap`
  with `~proplot.styletools.PerceptuallyUniformColormap` (:commit:`972956b1`).
* Include the `alpha` channel when saving colormaps and cycles by default
  (:pr:`149`, :commit:`117e05f2`).
* Permit 8-character hex strings with alpha channels when loading colormaps
  and color cycles from hex files (:pr:`149`, :commit:`381a84d4`).
* Publicly support "filling" axes with colorbars using ``loc='fill'``
  (:pr:`149`, :commit:`057c9895`).
* Return both figure and axes in ``show_`` functions; this gives users access
  to the axes and prevents drawing them twice in notebooks
  (:pr:`149`, :commit:`2f600bc9`).
* Enable passing callables to `~proplot.axistools.Formatter` to create a
  `~proplot.axistools.FuncFormatter` instance.
* Support sampling `~prolot.styletools.LinearSegmentedColormap` into
  `~proplot.styletools.ListedColormaps` inside of
  `~proplot.styletools.Colormap` rather than `~proplot.styletools.Cycle`
  (:issue:`84`, :commit:`972956b1`).

Bug fixes
---------

* Fix various issues with axis label sharing and axis sharing for
  twinned axes and panel axes (:pr:`164`).
* Permit modifying existing cartopy geographic features with successive
  calls to `~proplot.axes.GeoAxes.format` (:pr:`168`).
* Fix issue drawing bar plots with datetime *x* axes (:pr:`156`).
* Fix issue where `~proplot.ticker.AutoFormatter` tools were not locale-aware, i.e. use
  comma as decimal point sometimes (:pr:`152`, :commit:`c7636296`).
* Fix issue where `~proplot.ticker.AutoFormatter` nonzero-value correction algorithm was
  right for wrong reasons and could be wrong in rare circumstances
  (:pr:`152`, :commit:`c7636296`).
* Fix issue where ``matplotlib.style.use`` resets backend
  (:pr:`149`, :commit:`c8319104`).
* Fix issue with colormaps with dots in name (:pr:`149`, :commit:`972956b1`).
* Fix logarithmic scale argument parsing deprecation (:pr:`149`, :commit:`6ed7dbc5`).
* Fix deprecation of direct access to ``matplotlib.cm.cmap_d``
  in matplotlib >=3.2.0 (:pr:`149`, :commit:`a69c16da`).
* Fix issues with string font sizes (:pr:`149`, :commit:`6121de03`). Add hidden
  `~proplot.config.RcConfigurator._get_font_size` method to
  translate font size to numeric.
* Fix issue where passing actual projection instances generated with
  `~proplot.constructor.Proj` to `~proplot.ui.subplots` could incorrectly
  pair cartopy projections with basemap axes and vice versa (:pr:`149`).
* Fix issue where could not draw colorbar from list of single-color
  `~matplotlib.collections.PathCollection`\ s, i.e.
  scatter plots (:pr:`149`, :commit:`e893900b`).
* Fix issue where importing proplot in jupyter notebooks resets the default
  inline backend (:pr:`149`, :commit:`6121de03`).
* Improve axis label sharing algorithm (:commit:`6535b219`).
* Fix main axis label sharing bugs in presence of panels
  (:commit:`7b709db9`).
* Fix v0.4.0 regression where panel sharing no longer works
  (:commit:`289e5538`).
* Fix `~proplot.axistools.AutoFormatter` bug with values close
  to zero (:issue:`124`, :commit:`9b7f89fd`)
* Fix `~proplot.axistools.AutoFormatter` bug with small negative
  numbers (:issue:`117`).
* Label cyclic Scientific colour maps as cyclic (:commit:`e10a3109`).
* Permit special colormap normalization and level scaling for
  colormap-colored contour plots, just like contourf (:pr:`149`, :commit:`054cceb5`).

Internals
---------

* **Major** internal change: Move functions into smaller separate
  files to mimic how matplotlib library is divided up (:pr:`149`).
* Add `internals` folder containing default proplot rc params, deprecation
  helper functions, and other internal tools (:pr:`149`).
* Make colorbar axes instances of `~proplot.axes.CartesianAxes`, just
  like panel axes.
* Rename ubiquitous `_notNone` function to `_not_none` and change to more
  sensible behavior.
* Turn some private `~proplot.config` functions into static
  methods (:commit:`6121de03`).
* Remove "smart bounds" feature from `FuncScale` (:pr:`166`, :commit:`9ac149ea`).
* Clean up axes iterators (:pr:`149`, :commit:`c8a0768a`).

Documentation
-------------

* Call figure objects `fig` instead of `f`.
* Major clean up of notebook examples (:commit:`f86542b5`).
* Major clean up `~proplot.wrappers` documentation (:commit:`9648c18f`)
* Fix dead "See Also" links (:commit:`d32c6506`).
* Use "Other parameters" tables more often (:commit:`d32c6506`).


ProPlot v0.5.0 (2020-02-10)
===========================

Deprecated
----------

* Remove `abcformat` from `~proplot.axes.Axes.format` (:commit:`2f295e18`).
* Rename `top` to `abovetop` in `~proplot.axes.Axes.format` (:commit:`500dd381`).
* Rename `abc.linewidth` and `title.linewidth` to ``borderwidth`` (:commit:`54eb4bee`).
* Rename `~proplot.wrappers.text_wrapper` `linewidth` and `invert` to
  `borderwidth` and `borderinvert` (:commit:`54eb4bee`).

Features
--------

* Add back `Fabio Crameri's scientific colour maps
  <http://www.fabiocrameri.ch/colourmaps.php>`__ (:pr:`116`).
* Permit both e.g. `locator` and `xlocator` as keyword arguments to
  `~proplot.axes.Axes.altx`, etc. (:commit:`57fab860`).
* Permit *descending* `~proplot.styletools.BinNorm` and
  `~proplot.styletools.LinearSegmentedNorm` levels (:pr:`119`).
* Permit overriding the font weight, style, and stretch in the
  `~proplot.styletools.show_fonts` table (:commit:`e8b9ee38`).
* Permit hiding "unknown" colormaps and color cycles in the
  `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  tables (:commit:`cb206f19`).

Bug fixes
---------

* Fix issue where `~proplot.styletools.show_cmaps` and
  `~proplot.styletools.show_cycles` colormap names were messed up
  (:commit:`13045599`)
* Fix issue where `~proplot.styletools.show_cmaps` and
  `~proplot.styletools.show_cycles` did not return figure instance
  (:commit:`98209e87`).
* Fix issue where user `values` passed to
  `~proplot.wrappers.colorbar_wrapper` were sometimes ignored
  (:commit:`fd4f8d5f`).
* Permit passing *lists of colors* to manually shade line contours and filled
  contours in `~proplot.wrappers.cmap_changer`.
* Prevent formatting rightmost meridian label as ``1e-10`` on cartopy map
  projections (:commit:`37fdd1eb`).
* Support CF-time axes by fixing bug in `~proplot.wrappers.standardize_1d`
  and `~proplot.wrappers.standardize_2d` (:issue:`103`, :pr:`121`).
* Redirect to the "default" location when using ``legend=True`` and
  ``colorbar=True`` to generate on-the-fly legends and colorbars
  (:commit:`c2c5c58d`). This feature was accidentally removed.
* Let `~proplot.wrappers.colorbar_wrapper` accept lists of colors
  (:commit:`e5f11591`). This feature was accidentally removed.

Internals
---------

* Remove various unused keyword arguments (:commit:`33654a42`).
* Major improvements to the API controlling axes titles and a-b-c labels
  (:commit:`1ef7e65e`).
* Always use full names ``left``, ``right``, ``top``, and ``bottom`` instead
  of ``l``, ``r``, ``b``, and ``t``, for clarity (:commit:`1ef7e65e`).
* Improve ``GrayCycle`` colormap, is now much shorter and built from
  reflected Fabio ``GrayC`` colormaps (:commit:`5b2c7eb7`).


ProPlot v0.4.3 (2020-01-21)
===========================

Deprecated
----------

* Remove `~proplot.rctools.ipython_autoreload`,
  `~proplot.rctools.ipython_autosave`, and `~proplot.rctools.ipython_matplotlib`
  (:issue:`112`, :pr:`113`). Move inline backend configuration to a hidden
  method that gets called whenever the ``rc_configurator`` is initalized.

Features
--------

* Permit comments at the head of colormap and color files
  (:commit:`0ffc1d15`).
* Make `~proplot.axes.Axes.parametric` match ``plot`` autoscaling behavior
  (:commit:`ecdcba82`).

Internals
---------

* Use `~proplot.axes.Axes.colorbar` instead of `~matplotlib.axes.Axes.imshow`
  for `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  displays (:pr:`107`).

ProPlot v0.4.2 (2020-01-09)
===========================

Features
--------

* Add ``family`` keyword arg to `~proplot.styletools.show_fonts` (:pr:`106`).
* Package the `TeX Gyre <http://www.gust.org.pl/projects/e-foundry/tex-gyre>`__
  font series with ProPlot (:pr:`106`). Remove a couple other fonts.
* Put the TeX Gyre fonts at the head of the serif, sans-serif, monospace,
  cursive, and fantasy ``rcParams`` font family lists (:issue:`104`, :pr:`106`).

Bug fixes
---------

* Fix issues with Fira Math weights unrecognized by matplotlib (:pr:`106`).

ProPlot v0.4.1 (2020-01-08)
===========================

Deprecation
-----------

* Change the default ``.proplotrc`` format from YAML to the ``.matplotlibrc``
  syntax (:pr:`101`).

Features
--------

* Comments (lines starting with ``#``) are now permitted in all RGB and HEX style
  colormap and cycle files (:pr:`100`).
* Break down `~proplot.styletools.show_cycles` bars into categories, just
  like `~proplot.styletools.show_cmaps` (:pr:`100`).

Bug fixes
---------

* Fix issue where `~proplot.styletools.show_cmaps` and `~proplot.styletools.show_cycles`
  draw empty axes (:pr:`100`).
* Add back the :ref:`default .proplorc file <The .proplotrc file>` to docs (:pr:`101`).
  To do this, ``conf.py`` auto-generates a file in ``_static``.

Internals
---------

* Add ``geogrid.color/linewidth/etc`` and ``gridminor.color/linewidth/etc``
  props as *children* of ``grid.color/linewidth/etc`` (:pr:`101`).
* Various `~proplot.rctools.rc_configurator` improvements, remove outdated
  global variables (:pr:`101`).
* Better error handling when loading colormap/cycle files, and calls to
  `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` now raise
  errors while calls to `~proplot.styletools.register_cmaps` and
  `~proplot.styletools.register_cycles` still issue warnings (:pr:`100`).

ProPlot v0.4.0 (2020-01-07)
===========================

Deprecated
----------

* Rename `basemap_defaults` to `~proplot.projs.basemap_kwargs` and
  `cartopy_projs` to `~proplot.projs.cartopy_names` (:commit:`431a06ce`).
* Remove ``subplots.innerspace``, ``subplots.titlespace``,
  ``subplots.xlabspace``, and ``subplots.ylabspace`` spacing arguments,
  automatically calculate default non-tight spacing using `~proplot.subplots._get_space`
  based on current tick lengths, label sizes, etc.
* Remove redundant `~proplot.rctools.use_fonts`, use
  ``rcParams['sans-serif']`` precedence instead (:pr:`95`).
* `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy` no longer accept
  "scale-spec" arguments.  Must be a function, two functions, or an axis
  scale instance (:pr:`96`).
* Remove `~proplot.axes.Axes` ``share[x|y]``, ``span[x|y]``, and
  ``align[x|y]`` kwargs (:pr:`99`).  These settings are now always
  figure-wide.
* Rename `~proplot.styletools.Cycle` ``samples`` to ``N``, rename
  `~proplot.styletools.show_colors` ``nbreak`` to ``nhues`` (:pr:`98`).

Features
--------

* Add `~proplot.styletools.LinearSegmentedColormap.from_file` static methods
  (:pr:`98`).  You can now load files by passing a name to
  `~proplot.styletools.Colormap`.
* Add TeX Gyre Heros as open source Helvetica-alternative; this is the new
  default font.  Add Fira Math as DejaVu Sans-alternative; has complete set
  of math characters (:pr:`95`).
* Add `xlinewidth`, `ylinewidth`, `xgridcolor`, `ygridcolor` keyword args to
  `~proplot.axes.XYAxes.format` (:pr:`95`).
* Add getters and setters for various `~proplot.subplots.Figure` settings
  like ``share[x|y]``, ``span[x|y]``, and ``align[x|y]`` (:pr:`99`).
* Let `~proplot.axes.Axes.twinx`, `~proplot.axes.Axes.twiny`,
  `~proplot.axes.Axes.altx`, and `~proplot.axes.Axes.alty` accept
  `~proplot.axes.XYAxes.format` keyword args just like
  `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy` (:pr:`99`).
* Add `~proplot.subplots.Figure` ``fallback_to_cm`` kwarg. This is used by
  `~proplot.styletools.show_fonts` to show dummy glyphs to clearly illustrate
  when fonts are missing characters, but preserve graceful fallback for end
  user.
* Improve `~proplot.projs.Proj` constructor function. It now accepts
  `~cartopy.crs.Projection` and `~mpl_toolkits.basemap.Basemap` instances,
  just like other constructor functions, and returns only the projection
  instance (:pr:`92`).
* `~proplot.rctools.rc` `~proplot.rctools.rc_configurator.__getitem__` always
  returns the setting. To get context block-restricted settings, you must
  explicitly pass ``context=True`` to `~proplot.rctools.rc_configurator.get`,
  `~proplot.rctools.rc_configurator.fill`, or
  `~proplot.rctools.rc_configurator.category` (:pr:`91`).

Bug fixes
---------

* Fix `~proplot.rctools.rc_configurator.context` bug (:issue:`80` and :pr:`91`).
* Fix issues with `~proplot.axes.Axes.dualx` and `~proplot.axes.Axes.dualy`
  with non-linear parent scales (:pr:`96`).
* Ignore TTC fonts because they cannot be saved in EPS/PDF figures
  (:issue:`94` and :pr:`95`).
* Do not try to use Helvetica Neue because "thin" font style is read as
  regular (:issue:`94` and :pr:`95`).

Documentation
-------------

* Use the imperative mood for docstring summaries (:pr:`92`).
* Fix `~proplot.styletools.show_cycles` bug (:pr:`90`) and show cycles using
  colorbars rather than lines (:pr:`98`).

Internals
---------

* Define `~proplot.rctools.rc` default values with inline dictionaries rather
  than with a default ``.proplotrc`` file, change the auto-generated user
  ``.proplotrc`` (:pr:`91`).
* Remove useless `panel_kw` keyword arg from
  `~proplot.wrappers.legend_wrapper` and `~proplot.wrappers.colorbar_wrapper`
  (:pr:`91`). Remove `wflush`, `hflush`, and `flush` keyword args from
  `~proplot.subplots.subplots` that should have been removed long ago.

ProPlot v0.3.1 (2019-12-16)
===========================

Bug fixes
---------

* Fix issue where custom fonts were not synced (:commit:`a1b47b4c`).
* Fix issue with latest versions of matplotlib where ``%matplotlib inline``
  fails *silently* so the backend is not instantiated (:commit:`cc39dc56`).

ProPlot v0.3.0 (2019-12-15)
===========================

Deprecated
----------

* Remove ``'Moisture'`` colormap (:commit:`cf8952b1`).

Features
--------

* Add `~proplot.styletools.use_font`, only sync Google Fonts fonts
  (:pr:`87`).
* New ``'DryWet'`` colormap is colorblind friendly (:commit:`0280e266`).
* Permit shifting arbitrary colormaps by ``180`` degrees by appending the
  name with ``'_shifted'``, just like ``'_r'`` (:commit:`e2e2b2c7`).

Bug fixes
---------

* Add brute force workaround for saving colormaps with *callable* segmentdata
  (:commit:`8201a806`).
* Fix issue with latest versions of matplotlib where ``%matplotlib inline``
  fails *silently* so the backend is not instantiated (:commit:`cc39dc56`).
* Fix `~proplot.styletools.LinearSegmentedColormap.shifted` when `shift` is
  not ``180`` (:commit:`e2e2b2c7`).
* Save the ``cyclic`` and ``gamma`` attributes in JSON files too
  (:commit:`8201a806`).

Documentation
-------------

* Cleanup notebooks, especially the colormaps demo (e.g. :commit:`952d4cb3`).

Internals
---------

* Change `~time.clock` to `~time.perf_counter` (:pr:`86`).

ProPlot v0.2.7 (2019-12-09)
===========================

Bug fixes
---------

* Fix issue where `~proplot.styletools.AutoFormatter` logarithmic scale
  points are incorrect (:commit:`9b164733`).

Documentation
-------------

* Improve :ref:`Configuring proplot` documentation (:commit:`9d50719b`).

Internals
---------

* Remove `prefix`, `suffix`, and `negpos` keyword args from
  `~proplot.styletools.SimpleFormatter`, remove `precision` keyword arg from
  `~proplot.styletools.AutoFormatter` (:commit:`8520e363`).
* Make ``'deglat'``, ``'deglon'``, ``'lat'``, ``'lon'``, and ``'deg'``
  instances of `~proplot.styletools.AutoFormatter` instead of
  `~proplot.styletools.SimpleFormatter` (:commit:`8520e363`). The latter
  should just be used for contours.

ProPlot v0.2.6 (2019-12-08)
===========================

Bug fixes
---------

* Fix issue where twin axes are drawn *twice* (:commit:`56145122`).


ProPlot v0.2.5 (2019-12-07)
===========================

Features
--------

* Much better `~proplot.axistools.CutoffScale` algorithm, permit arbitrary
  cutoffs (:pr:`83`).

ProPlot v0.2.4 (2019-12-07)
===========================

Deprecated
----------

* Rename `ColorCacheDict` to `~proplot.styletools.ColorDict`
  (:commit:`aee7d1be`).
* Rename `colors` to `~proplot.styletools.Colors` (:commit:`aee7d1be`)
* Remove `fonts_system` and `fonts_proplot`, rename `colordict` to
  `~proplot.styletools.colors`, make top-level variables more robust
  (:commit:`861583f8`).

Documentation
-------------

* Params table for `~proplot.styletools.show_fonts` (:commit:`861583f8`).

Internals
---------

* Improvements to `~proplot.styletools.register_colors`.

ProPlot v0.2.3 (2019-12-05)
===========================

Bug fixes
---------

* Fix issue with overlapping gridlines using monkey patches on gridliner
  instances (:commit:`8960ebdc`).
* Fix issue where auto colorbar labels are not applied when ``globe=True``
  (:commit:`ecb3c899`).
* More sensible zorder for gridlines (:commit:`90d94e55`).
* Fix issue where customized super title settings are overridden when new
  axes are created (:commit:`35cb21f2`).

Documentation
-------------

* Organize ipython notebook documentation (:commit:`35cb21f2`).

Internals
---------

* Major cleanup of the `~proplot.wrappers.colorbar_wrapper` source code,
  handle minor ticks using the builtin matplotlib API just like major ticks
  (:commit:`b9976220`).

ProPlot v0.2.2 (2019-12-04)
===========================

Deprecated
----------

* Rename `~proplot.subplots.axes_grid` to `~proplot.subplots.subplot_grid`
  (:commit:`ac14e9dd`).

Bug fixes
---------

* Fix shared *x* and *y* axis bugs (:commit:`ac14e9dd`).

Documentation
-------------

* Make notebook examples PEP8 compliant (:commit:`97f5ffd4`). Much more
  readable now.

ProPlot v0.2.1 (2019-12-02)
===========================

Deprecated
----------

* Rename `autoreload_setup`, `autosave_setup`, and `matplotlib_setup` to
  `~proplot.rctools.ipython_autoreload`, `~proplot.rctools.ipython_autosave`,
  and `~proplot.rctools.ipython_matplotlib`, respectively
  (:commit:`84e80c1e`).

ProPlot v0.2.0 (2019-12-02)
===========================

Deprecated
----------

* Remove the ``nbsetup`` rc setting in favor of separate ``autosave``,
  ``autoreload``, and ``matplotlib`` settings for triggering the respective
  ``%`` magic commands.  (:commit:`3a622887`; ``nbsetup`` is still accepted
  but no longer documented).
* Rename the ``format`` rc setting in favor of the ``inlinefmt`` setting
  (:commit:`3a622887`; ``format`` is still accepted but no longer
  documented).
* Rename ``FlexibleGridSpec`` and ``FlexibleSubplotSpec`` to ``GridSpec`` and
  ``SubplotSpec`` (:commit:`3a622887`; until :pr:`110` is merged it is
  impossible to use these manually, so this won't bother anyone).

Features
--------

* Support manual resizing for all backends, including ``osx`` and ``qt``
  (:commit:`3a622887`).

Bug fixes
---------

* Disable automatic resizing for the ``nbAgg`` interactive inline backend.
  Found no suitable workaround (:commit:`3a622887`).

Internals
---------

* Organize the ``rc`` documentation and the default ``.proplotrc`` file
  (:commit:`3a622887`).
* Rename ``rcParamsCustom`` to ``rcParamsLong`` (:commit:`3a622887`; this is
  inaccessible to the user).
* When calling ``fig.canvas.print_figure()`` on a stale figure, call
  ``fig.canvas.draw()`` first. May be overkill for
  `~matplotlib.figure.Figure.savefig` but critical for correctly displaying
  already-drawn notebook figures.

ProPlot v0.1.0 (2019-12-01)
===========================

Internals
---------

* Include `flake8` in Travis CI testing (:commit:`8743b857`).
* Enforce source code PEP8 compliance (:commit:`78da51a7`).
* Use pre-commit for all future commits (:commit:`e14f6809`).
* Implement tight layout stuff with canvas monkey patches
  (:commit:`67221d10`).  ProPlot now works for arbitrary backends, not just
  inline and qt.

Documentation
-------------

* Various `RTD bugfixes
  <https://github.com/readthedocs/readthedocs.org/issues/6412>`__ (e.g.
  :commit:`37633a4c`).

ProPlot v0.0.0 (2019-11-27)
===========================

The first version released on `PyPi <https://pypi.org/project/proplot/>`__.

.. _Luke Davis: https://github.com/lukelbd

.. _Riley Brady: https://github.com/bradyrx

.. _Stephane Raynaud: https://github.com/stefraynaud

.. _Mickaël Lalande: https://github.com/mickaellalande

.. _Pratiman Patel: https://github.com/pratiman-91

.. _Zachary Moon: https://github.com/zmoon

.. _Eli Knaap: https://github.com/knaaptime
