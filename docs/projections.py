# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
#
# .. _polar: https://matplotlib.org/3.1.0/gallery/pie_and_polar_charts/polar_demo.html
#
# .. _cartopy: https://scitools.org.uk/cartopy/docs/latest/
#
# .. _basemap: https://matplotlib.org/basemap/index.html
#
# .. _ug_proj:
#
# Geographic and polar plots
# ==========================
#
# ProPlot includes several advanced features for working with `polar`_
# and :ref:`geographic projections <ug_geo>`.
#
# To change the axes projection, pass ``proj='name'`` to an axes-creation command
# (i.e., `~proplot.figure.Figure.add_subplot`, `~proplot.figure.Figure.add_subplots`,
# `~proplot.figure.Figure.subplot`, or `~proplot.figure.Figure.subplots`). To use
# different projections for different subplots when creating your subplots al
# at once with `~proplot.figure.Figure.subplots`, pass either a list of projection
# names or a dictionary of projection names with the subplot number as the key.
# For example, a 2-column figure with a Cartesian axes on the left and a Plate Carrée
# projection on the right can be built with either ``proj=('cartesian', 'pcarree')``
# or ``proj={2: 'pcarree'}``. The default projection is `~proplot.axes.CartesianAxes`,
# optionally specified with the key ``'cartesian'``.

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geo:
#
# Geographic axes
# ---------------

# To create geographic axes, pass e.g. ``proj='name'`` to an
# axes-creation command (see :ref:`above <ug_proj>`) where ``name`` is any valid
# :ref:`PROJ projection name <proj_included>`. Alternatively, you can use
# ``proj=projection_instance`` where ``projection_instance`` is an object returned
# by the `~proplot.constructor.Proj` :ref:`constructor function <why_constructor>`
# (see below for details). Requesting geographic projections returns
# `~proplot.axes.GeoAxes` instance(s) with their own `~proplot.axes.GeoAxes.format`
# command. `proplot.axes.GeoAxes.format` facilitates :ref:`geographic-specific
# modifications <ug_geoformat>` like meridional and parallel gridlines and land
# mass outlines. The syntax is very similar to `proplot.axes.CartesianAxes.format`.
# In the below example, we create and format a very simple geographic plot.

# %%
# Use an on-the-fly projection
import proplot as pplt
fig = pplt.figure(refwidth=3)
axs = fig.subplots(nrows=2, proj='robin', proj_kw={'lon_0': 180})
# proj = pplt.Proj('robin', lon_0=180)
# axs = pplt.subplots(nrows=2, proj=proj)  # equivalent to above
axs.format(
    suptitle='Figure with single projection',
    land=True, latlines=30, lonlines=60,
)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_backends:
#
# Cartopy and basemap
# -------------------
#
# The `proplot.axes.GeoAxes` class uses either `cartopy`_ or `basemap`_ as "backends"
# to :ref:`format the axes <ug_geoformat>` and :ref:`plot stuff <ug_geoplot>` in
# the axes. A few details:
#
# * Cartopy is the default backend. When you request projection names with cartopy
#   as the backend (or pass a `cartopy.crs.Projection` to the `proj` keyword), the
#   returned axes is a subclass of `cartopy.mpl.geoaxes.GeoAxes`. Under the hood,
#   invoking `~proplot.axes.GeoAxes.format` with cartopy as the backend changes map
#   bounds using `~cartopy.mpl.geoaxes.GeoAxes.set_extent`, adds major and minor
#   gridlines using `~cartopy.mpl.geoaxes.GeoAxes.gridlines`, and adds geographic
#   features using `~cartopy.mpl.geoaxes.GeoAxes.add_feature`. If you prefer, you can
#   use the standard `cartopy.mpl.geoaxes.GeoAxes` methods just like you would in
#   cartopy. If you need to use the underlying `~cartopy.crs.Projection` instance, it
#   is available via the `~proplot.axes.GeoAxes.projection` attribute.
#
# * Basemap is an alternative backend. To use basemap, set :rcraw:`basemap` to
#   ``True`` or pass ``basemap=True`` to the axes-creation command. When you
#   request a projection name with basemap as the backend (or pass a
#   `~mpl_toolkits.basemap.Basemap` to the `proj` keyword), the returned axes
#   redirects the plotting methods plot, scatter, contour, contourf, pcolor,
#   pcolormesh, quiver, streamplot, and barb to the identically named methods on
#   the `~mpl_toolkits.basemap.Basemap` instance. This means you can work
#   with the standard axes plotting methods rather than the basemap methods --
#   just like cartopy. Under the hood, invoking `~proplot.axes.GeoAxes.format`
#   with basemap as the backend adds major and minor gridlines using
#   `~mpl_toolkits.basemap.Basemap.drawmeridians` and
#   `~mpl_toolkits.basemap.Basemap.drawparallels` and adds geographic features
#   using methods like `~mpl_toolkits.basemap.Basemap.fillcontinents`
#   and `~mpl_toolkits.basemap.Basemap.drawcoastlines`. If you need to
#   use the underlying `~mpl_toolkits.basemap.Basemap` instance, it is
#   available via the `~proplot.axes.GeoAxes.projection` attribute.
#
# Together, these features let you work with geophysical data without invoking
# verbose cartopy classes like `~cartopy.crs.LambertAzimuthalEqualArea` or
# keeping track of separate `~mpl_toolkits.basemap.Basemap` instances. This
# considerably reduces the amount of code needed to make complex geographic
# plots. In the below examples, we create a variety of plots using both
# cartopy and basemap as backends.
#
# .. note::
#
#    * By default, ProPlot gives circular boundaries to polar cartopy projections
#      like `~cartopy.crs.NorthPolarStereo` (see `this example
#      <https://scitools.org.uk/cartopy/docs/latest/gallery/lines_and_polygons/always_circular_stereo.html>`__
#      from the cartopy website). This is consistent with basemap's default behavior.
#      To disable this feature, set :rcraw:`cartopy.circular` to ``False``.
#      Please note that cartopy cannot add gridline labels to polar plots
#      with circular boundaries.
#    * By default, ProPlot uses `~cartopy.mpl.geoaxes.GeoAxes.set_global` to give
#      non-polar cartopy projections global extent and bounds polar cartopy projections
#      at the equator. This is a deviation from cartopy, which determines map boundaries
#      automatically based on the coordinates of the plotted content. To revert to
#      cartopy's default behavior, set :rcraw:`cartopy.autoextent` to ``True``.
#    * To make things more consistent, the `~proplot.constructor.Proj` constructor
#      function lets you supply native `PROJ <https://proj.org>`__ keyword names
#      for the cartopy `~cartopy.crs.Projection` classes (e.g., `lon_0` instead
#      of `central_longitude`) and instantiates `~mpl_toolkits.basemap.Basemap`
#      projections with sensible default PROJ parameters rather than raising an error
#      when they are omitted (e.g., ``lon_0=0`` as the default for most projections).
#    * Basemap is `no longer maintained \
#      <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__
#      and will not work with matplotlib versions more recent than 3.2.2. However,
#      basemap gridline labels often look nicer than cartopy -- especially when
#      "inline" cartopy labels are disabled. This is the main reason ProPlot continues
#      to support basemap. When cartopy gridline labels improve, basemap support
#      may be deprecated.

# %%
import proplot as pplt
fig = pplt.figure()

# Add projections
gs = pplt.GridSpec(ncols=2, nrows=3, hratios=(1, 1, 1.4))
for i, proj in enumerate(('cyl', 'hammer', 'npstere')):
    ax1 = fig.subplot(gs[i, 0], proj=proj, basemap=True)  # basemap
    ax2 = fig.subplot(gs[i, 1], proj=proj)  # cartopy

# Format projections
fig.format(
    land=True,
    suptitle='Figure with several projections',
    toplabels=('Basemap projections', 'Cartopy projections'),
    toplabelweight='normal',
    latlines=30, lonlines=60,
    lonlabels='b', latlabels='r',  # or lonlabels=True, labels=True, etc.
)
fig.subplotgrid[-2:].format(latlines=20, lonlines=30)  # dense gridlines for polar plots
pplt.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geoplot:
#
# Plotting in projections
# -----------------------
#
# In ProPlot, plotting with `~proplot.axes.GeoAxes` is very similar to plotting
# with `~proplot.axes.CartesianAxes`. ProPlot makes longitude-latitude
# (i.e., Plate Carrée) coordinates the *default* coordinate system by passing
# ``transform=ccrs.PlateCarree()`` to cartopy plotting commands and ``latlon=True`` to
# basemap plotting commands. And again, note that basemap plotting commands are invoked
# from the `proplot.axes.GeoAxes` rather than the `~mpl_toolkits.basemap.Basemap`
# instance -- just like cartopy. When using basemap as the "backend", you should not
# have to work with the `~mpl_toolkits.basemap.Basemap` instance directly.
#
# To ensure that graphics generated by :ref:`plotting commands <ug_2dplots>` like
# `~matplotlib.axes.Axes.contour` fill the entire globe, simply pass ``globe=True``
# to the command. This interpolates the data to the poles and across the longitude
# seam before plotting. This is a convenient and succinct alternative to cartopy's
# `~cartopy.util.add_cyclic_point` and basemap's `~mpl_toolkits.basemap.addcyclic`.
#
# Geographic features can be drawn underneath data or on top of data by changing the
# corresponding `zorder <https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
# setting. For example, to draw land patches on top of all plotted content as
# a "land mask," use ``ax.format(land=True, landzorder=4)`` or set :rcraw:`land.zorder`
# to ``True``. See the :ref:`next section <ug_geoformat>` for details.

# %%
import proplot as pplt
import numpy as np

# Fake data with unusual longitude seam location and without coverage over poles
offset = -40
lon = pplt.arange(offset, 360 + offset - 1, 60)
lat = pplt.arange(-60, 60 + 1, 30)
state = np.random.RandomState(51423)
data = state.rand(len(lat), len(lon))

# Plot data both without and with globe=True
for globe in (False, True):
    string = 'with' if globe else 'without'
    gs = pplt.GridSpec(nrows=2, ncols=2)
    fig = pplt.figure(refwidth=2.5)
    for i, ss in enumerate(gs):
        ax = fig.subplot(ss, proj='kav7', basemap=(i % 2))
        cmap = ('sunset', 'sunrise')[i % 2]
        if i > 1:
            ax.pcolor(lon, lat, data, cmap=cmap, globe=globe, extend='both')
        else:
            m = ax.contourf(lon, lat, data, cmap=cmap, globe=globe, extend='both')
            fig.colorbar(m, loc='b', span=i + 1, label='values', extendsize='1.7em')
    fig.format(
        suptitle=f'Geophysical data {string} global coverage',
        toplabels=('Cartopy example', 'Basemap example'),
        leftlabels=('Filled contours', 'Grid boxes'),
        toplabelweight='normal', leftlabelweight='normal',
        coast=True, lonlines=90,
        abc='A.', abcloc='ul', abcborder=False,
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geoformat:
#
# Formatting projections
# ----------------------
#
# The `proplot.axes.GeoAxes.format` command facilitates geographic-specific axes
# modifications. It can toggle and configure the "major" and "minor" longitude and
# latitude gridline locations using the `grid`, `lonlocator`, `latlocator`, `gridminor`,
# `lonminorlocator`, and `latminorlocator` keys, and configure gridline label formatting
# with `lonformatter` and `latformatter` (analogous to `xlocator`, `xminorlocator`,
# and `xformatter` used by `proplot.axes.CartesianAxes.format`). By default, inline
# cartopy labels and cartopy label rotation are turned off, but inline labels can
# be turned on using ``loninline=True``, ``latinline=True``, or ``inlinelabels=True``
# or by setting :rcraw:`grid.inlinelabels` to ``True``, and label rotation can be
# turned on using ``rotatelabels=True`` or by setting :rcraw:`grid.rotatelabels`
# to ``True``. The padding between the map edge and the labels can be changed
# using `labelpad` or by changing :rcraw:`grid.labelpad`.

# `proplot.axes.GeoAxes.format` can also set the cartopy projection bounding longitudes
# and latitudes with `lonlim` and `latlim` (analogous to `xlim` and `ylim`), set the
# latitude bound for circular polar projections using `boundinglat`, and toggle and
# configure geographic features like land masses, coastlines, and administrative
# borders using :ref:`settings <rc_proplot>` like `land`, `landcolor`, `coast`,
# `coastcolor`, and `coastlinewidth`. Finally, since `proplot.axes.GeoAxes.format`
# calls `proplot.axes.Axes.format`, it can be used to add axes titles, a-b-c labels,
# and figure titles, just like `proplot.axes.CartesianAxes.format`.
#
# For details, see the `proplot.axes.GeoAxes.format` documentation.

# %%
import proplot as pplt
gs = pplt.GridSpec(ncols=3, nrows=2, wratios=(1, 1, 1.2), hratios=(1, 1.2))
fig = pplt.figure(refwidth=4)

# Styling projections in different ways
ax = fig.subplot(gs[0, :2], proj='eqearth')
ax.format(
    title='Equal earth', land=True, landcolor='navy', facecolor='pale blue',
    coastcolor='gray5', borderscolor='gray5', innerborderscolor='gray5',
    gridlinewidth=1.5, gridcolor='gray5', gridalpha=0.5,
    gridminor=True, gridminorlinewidth=0.5,
    coast=True, borders=True, borderslinewidth=0.8,
)
ax = fig.subplot(gs[0, 2], proj='ortho')
ax.format(
    title='Orthographic', reso='med', land=True, coast=True, latlines=10, lonlines=15,
    landcolor='mushroom', suptitle='Projection axes formatting demo',
    facecolor='petrol', coastcolor='charcoal', coastlinewidth=0.8, gridlinewidth=1
)
ax = fig.subplot(gs[1, :], proj='wintri')
ax.format(
    land=True, facecolor='ocean blue', landcolor='bisque', title='Winkel tripel',
    lonlines=60, latlines=15,
    gridlinewidth=0.8, gridminor=True, gridminorlinestyle=':',
    lonlabels=True, latlabels='r', loninline=True,
    gridlabelcolor='gray8', gridlabelsize='med-large',
)
fig.format(
    suptitle='Projection axes formatting demo',
    toplabels=('Column 1', 'Column 2'),
    abc='A.', abcloc='ul', abcborder=False, linewidth=1.5
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_zoom:
#
# Zooming into projections
# ------------------------
#
# To zoom into cartopy projections, use
# `~cartopy.mpl.geoaxes.GeoAxes.set_extent` or pass `lonlim`,
# `latlim`, or `boundinglat` to `~proplot.axes.GeoAxes.format`. The `boundinglat`
# keyword controls the circular latitude boundary for North Polar and
# South Polar Stereographic, Azimuthal Equidistant, Lambert Azimuthal
# Equal-Area, and Gnomonic projections. By default, ProPlot tries to use the
# degree-minute-second cartopy locators and formatters made available in cartopy
# 0.18. You can switch from minute-second subintervals to traditional decimal
# subintervals by passing ``dms=False`` to `~proplot.axes.GeoAxes.format`
# or by setting :rcraw:`grid.dmslabels` to ``False``.
#
# To zoom into basemap projections, pass any of the `boundinglat`,
# `llcrnrlon`, `llcrnrlat`, `urcrnrlon`, `urcrnrlat`, `llcrnrx`, `llcrnry`,
# `urcrnrx`, `urcrnry`, `width`, or `height` keyword arguments to
# the `~proplot.constructor.Proj` constructor function either directly or via
# the `proj_kw` `~proplot.ui.subplots` keyword argument. You can also pass
# `lonlim` and `latlim` to `~proplot.constructor.Proj` and these arguments
# will be used for `llcrnrlon`, `llcrnrlat`, etc. You cannot zoom into basemap
# projections with `format` after they have already been created.

# %%
import proplot as pplt

# Plate Carrée map projection
pplt.rc.reso = 'med'  # use higher res for zoomed in geographic features
basemap = pplt.Proj('cyl', lonlim=(-20, 180), latlim=(-10, 50), basemap=True)
fig, axs = pplt.subplots(nrows=2, refwidth=5, proj=('cyl', basemap))
axs.format(
    land=True, labels=True, lonlines=20, latlines=20,
    gridminor=True, suptitle='Zooming into projections'
)
axs[0].format(lonlim=(-140, 60), latlim=(-10, 50), labels=True)
axs[0].format(title='Cartopy example')
axs[1].format(title='Basemap example')

# %%
import proplot as pplt

# Pole-centered map projections
basemap = pplt.Proj('npaeqd', boundinglat=60, basemap=True)
fig, axs = pplt.subplots(ncols=2, refwidth=2.7, proj=('splaea', basemap))
fig.format(suptitle='Zooming into polar projections')
axs.format(land=True, latmax=80)  # no gridlines poleward of 80 degrees
axs[0].format(boundinglat=-60, title='Cartopy example')
axs[1].format(title='Basemap example')

# %%
import proplot as pplt

# Zooming in on continents
fig = pplt.figure(refwidth=3)
ax = fig.subplot(121, proj='lcc', proj_kw={'lon_0': 0})
ax.format(lonlim=(-20, 50), latlim=(30, 70), title='Cartopy example')
proj = pplt.Proj('lcc', lon_0=-100, lat_0=45, width=8e6, height=8e6, basemap=True)
ax = fig.subplot(122, proj=proj)
ax.format(lonlines=20, title='Basemap example')
fig.format(suptitle='Zooming into specific regions', land=True)


# %%
import proplot as pplt

# Zooming in with cartopy degree-minute-second labels
pplt.rc.reso = 'hi'
fig = pplt.figure(refwidth=2.5)
ax = fig.subplot(121, proj='cyl')
ax.format(lonlim=(-7.5, 2), latlim=(49.5, 59))
ax = fig.subplot(122, proj='cyl')
ax.format(lonlim=(-6, -2), latlim=(54.5, 58.5))
fig.format(
    land=True, labels=True,
    borders=True, borderscolor='white',
    suptitle='Cartopy degree-minute-second labels',
)
pplt.rc.reset()

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _proj_included:
#
# Included projections
# --------------------
#
# The available `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__
# and `basemap <https://matplotlib.org/basemap/index.html>`__ projections are
# plotted below. See `~proplot.constructor.Proj` for a table of projection
# names with links to the relevant `PROJ <https://proj.org>`__ documentation.
#
# ProPlot uses the cartopy API to add the Aitoff, Hammer, Winkel Tripel, and
# Kavrisky VII projections (i.e., ``'aitoff'``, ``'hammer'``, ``'wintri'``,
# and ``'kav7'``), as well as North and South polar versions of the Azimuthal
# Equidistant, Lambert Azimuthal Equal-Area, and Gnomic projections (i.e.,
# ``'npaeqd'``, ``'spaeqd'``, ``'nplaea'``, ``'splaea'``, ``'npgnom'``, and
# ``'spgnom'``), modeled after the existing `~cartopy.crs.NorthPolarStereo`
# and `~cartopy.crs.SouthPolarStereo` projections.

# %%
import proplot as pplt

# Table of cartopy projections
projs = [
    'cyl', 'merc', 'mill', 'lcyl', 'tmerc',
    'robin', 'hammer', 'moll', 'kav7', 'aitoff', 'wintri', 'sinu',
    'geos', 'ortho', 'nsper', 'aea', 'eqdc', 'lcc', 'gnom',
    'npstere', 'nplaea', 'npaeqd', 'npgnom', 'igh',
    'eck1', 'eck2', 'eck3', 'eck4', 'eck5', 'eck6'
]
fig, axs = pplt.subplots(ncols=3, nrows=10, figwidth=7, proj=projs)
axs.format(
    land=True, reso='lo', labels=False,
    suptitle='Table of cartopy projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)

# %%
import proplot as pplt

# Table of basemap projections
projs = [
    'cyl', 'merc', 'mill', 'cea', 'gall', 'sinu',
    'eck4', 'robin', 'moll', 'kav7', 'hammer', 'mbtfpq',
    'geos', 'ortho', 'nsper',
    'vandg', 'aea', 'eqdc', 'gnom', 'cass', 'lcc',
    'npstere', 'npaeqd', 'nplaea'
]
fig, axs = pplt.subplots(ncols=3, nrows=8, basemap=True, figwidth=7, proj=projs)
axs.format(
    land=True, labels=False,
    suptitle='Table of basemap projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_polar:
#
# Polar axes
# ----------
#
# To create `polar axes <polar_>`_, pass ``proj='polar'`` to an axes-creation
# command (see :ref:`above <ug_proj>`). This returns `proplot.axes.PolarAxes`
# instance(s) with their own `~proplot.axes.PolarAxes.format` command.
# `proplot.axes.PolarAxes.format` facilitates polar-specific axes modifications
# like changing the central radius `r0`, the zero azimuth location `theta0`,
# and the positive azimuthal direction `thetadir`. It also supports toggling and
# configuring the "major" and "minor" gridline locations with `grid`, `rlocator`,
# `thetalocator`, `gridminor`, `rminorlocator`, and `thetaminorlocator` and formatting
# the gridline labels with `rformatter` and `thetaformatter` (analogous to `xlocator`,
# `xformatter`, and` `xminorlocator` used by `proplot.axes.CartesianAxes.format`),
# and creating "annular" or "sector" plots by changing the radial or azimuthal
# bounds `rlim` and `thetalim`. Finally, since `proplot.axes.PolarAxes.format`
# calls `proplot.axes.Axes.format`, it can be used to add axes titles, a-b-c
# labels, and figure titles.
#
# For details, see `proplot.axes.PolarAxes.format`.

# %%
import proplot as pplt
import numpy as np
N = 200
state = np.random.RandomState(51423)
x = np.linspace(0, 2 * np.pi, N)[:, None] + np.arange(5) * 2 * np.pi / 5
y = 100 * (state.rand(N, 5) - 0.3).cumsum(axis=0) / N
fig, axs = pplt.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], proj='polar')
axs.format(
    suptitle='Polar axes demo', linewidth=1, titlepad='1em',
    ticklabelsize=9, rlines=0.5, rlim=(0, 19),
)
for ax in axs:
    ax.plot(x, y, cycle='FlatUI', zorder=0, lw=3)

# Standard polar plot
axs[0].format(
    title='Normal plot', thetaformatter='tau',
    rlabelpos=225, rlines=pplt.arange(5, 30, 5),
    edgecolor='red8', tickpad='1em',
)

# Sector plot
axs[1].format(
    title='Sector plot', thetadir=-1, thetalines=90, thetalim=(0, 270), theta0='N',
    rlim=(0, 22), rlines=pplt.arange(5, 30, 5),
)

# Annular plot
axs[2].format(
    title='Annular plot', thetadir=-1, thetalines=20, gridcolor='red',
    r0=-20, rlim=(0, 22), rformatter='null', rlocator=2
)
