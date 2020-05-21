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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
#
# .. _polar: https://matplotlib.org/3.1.0/gallery/pie_and_polar_charts/polar_demo.html

# .. _cartopy: https://scitools.org.uk/cartopy/docs/latest/

# .. _basemap: https://matplotlib.org/basemap/index.html
#
# .. _ug_proj:
#
# Geographic and polar plots
# ==========================
#
# ProPlot includes features for working with `polar axes <polar_>`_
# and the `cartopy`_ and `basemap`_ map projection packages. These features
# are optional -- installation of cartopy and basemap are not required.
#
# To change the axes projection, pass ``proj='name'`` to
# `~proplot.ui.subplots`. To use different projections for different
# subplots, pass a dictionary of projection names with the subplot number as
# the key -- for example, ``proj={1: 'name'}``. The default "projection" is
# always `~proplot.axes.CartesianAxes`.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_polar:
#
# Polar axes
# ----------
#
# To draw `polar axes <polar_>`_, pass ``proj='polar'`` or e.g. ``proj={1:'polar'}``
# to `~proplot.ui.subplots`. This generates a `~proplot.axes.PolarAxes`
# instance with its own `proplot.axes.PolarAxes.format` command. This
# command permits polar-specific modifications like changing the central radius `r0`,
# the zero azimuth location `theta0`, and the positive azimuthal direction `thetadir`.
# It also supports changing the radial and azimuthal limits `rlim` and `thetalim`,
# which can be used to make sector plots and annular plots.
#
# For details, see `proplot.axes.PolarAxes.format`.

# %%
import proplot as plot
import numpy as np
N = 200
state = np.random.RandomState(51423)
x = np.linspace(0, 2 * np.pi, N)
y = 100 * (state.rand(N, 5) - 0.3).cumsum(axis=0) / N
plot.rc['axes.titlepad'] = '1em'  # default matplotlib offset is incorrect
fig, axs = plot.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], proj='polar')
axs.format(
    suptitle='Polar axes demo', linewidth=1,
    ticklabelsize=9, rlines=0.5, rlim=(0, 19),
)
for i in range(5):
    xi = x + i * 2 * np.pi / 5
    axs.plot(xi, y[:, i], cycle='FlatUI', zorder=0, lw=3)

# Standard polar plot
axs[0].format(
    title='Normal plot', thetaformatter='pi', rlines=5,
    rlabelpos=180, color='gray8', tickpad='1em'
)

# Sector plot
axs[1].format(
    title='Sector plot', thetadir=-1, thetalines=90, thetalim=(0, 270), theta0='N',
    rlim=(0, 22), rlines=5
)

# Annular plot
axs[2].format(
    title='Annular plot', thetadir=-1, thetalines=20, gridcolor='red',
    r0=0, rlim=(10, 22), rformatter='null', rlocator=2
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geo:
#
# Geographic axes
# ---------------
#
# ProPlot can turn any subplot into a geographic projection using
# the `cartopy`_ or `basemap`_ packages as "backends." The
# `~proplot.axes.GeoAxes` class and the `~proplot.constructor.Proj`
# constructor function ensure that ProPlot's syntax with cartopy as the
# "backend" is exactly the same as when basemap is the "backend".
# Cartopy is the default backend, but you can switch
# to basemap using ``basemap=True`` (see below).
#
# To turn a subplot into a geographic projection, pass
# ``proj='name'`` or e.g. ``proj={2: 'name'}``
# (:ref:`see above <ug_proj>`) to  `~proplot.ui.subplots`
# where ``name`` is any valid :ref:`PROJ projection name <proj_included>`.
# You can also generate a `cartopy.crs.Projection` or `mpl_toolkits.basemap.Basemap`
# instance directly using the `~proplot.constructor.Proj` constructor function and
# pass the class instance to `proj`.
#
# When you request map projections, `~proplot.ui.subplots` returns
# instances of `~proplot.axes.CartopyAxes` or `~proplot.axes.BasemapAxes`.
# These axes have the following properties:
#
# * `~proplot.axes.CartopyAxes` joins the cartopy
#   `~cartopy.mpl.geoaxes.GeoAxes` class with the ProPlot
#   `~matplotlib.axes.Axes` class and adds a `~proplot.axes.GeoAxes.format`
#   command. This class includes all the normal `~cartopy.mpl.geoaxes.GeoAxes`
#   methods, and its `~proplot.axes.GeoAxes.format` method can be used to set
#   the map bounds with `~cartopy.mpl.geoaxes.GeoAxes.set_extent` and add
#   geographic features with `~cartopy.mpl.geoaxes.GeoAxes.add_feature`.
#
# * `~proplot.axes.BasemapAxes` redirects the plot, scatter, contour,
#   contourf, pcolor, pcolormesh, quiver, streamplot, and barb *axes methods* to
#   identically named methods on the `~mpl_toolkits.basemap.Basemap` instance.
#   This means you can work with axes plotting methods just like in cartopy.
#   `~proplot.axes.BasemapAxes` also provides access to `~mpl_toolkits.basemap.Basemap`
#   geographic plotting commands like `~mpl_toolkits.basemap.Basemap.fillcontinents`
#   via `~proplot.axes.GeoAxes.format`.
#
# These features help address the limitations of the cartopy and basemap APIs.
# You no longer have to invoke verbose cartopy classes like
# `~cartopy.crs.LambertAzimuthalEqualArea` and `~cartopy.feature.NaturalEarthFeature`,
# and you no longer have to directly work with the `~mpl_toolkits.basemap.Basemap`
# instance. However if you do need access to the projection class instances,
# they are stored as `proplot.axes.CartopyAxes.projection` and
# `proplot.axes.BasemapAxes.projection` attributes. Also, to make things more
# consistent, the `~proplot.constructor.Proj` constructor function lets you supply
# native `PROJ <https://proj.org>`__ keyword names to the cartopy
# `~cartopy.crs.Projection` classes (e.g. `lon_0` instead of `central_longitude`),
# and instantiates `~mpl_toolkits.basemap.Basemap` projections with sensible
# defaults (e.g. ``lon_0=0``) rather than raising an error when projection
# arguments are omitted.
#
# .. note::
#
#   ProPlot makes sure polar cartopy projections like `~cartopy.crs.NorthPolarStereo`
#   have a circular boundary. By default, polar projections are bounded at the
#   equator and non-polar projections are forced to have global extent with
#   `~cartopy.mpl.geoaxes.GeoAxes.set_global`. To revert to the behavior where
#   cartopy automatically determines map boundaries based on plotted content,
#   simply set :rcraw:`cartopy.autoextent` to ``True`` or
#   pass ``autoextent=True`` to `~proplot.axes.CartopyAxes`. See the
#   :ref:`configuration guide <ug_config>` for details.
#
# .. warning::
#
#    Basemap is `no longer a maintained package\
#    <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__.
#    However as shown below, gridline labels tend to look much nicer in basemap
#    than in cartopy -- especially when "inline" cartopy labels are disabled.
#    This is the main reason ProPlot continues to support both basemap and cartopy.
#    When cartopy catches up, basemap support may be deprecated.

# %%
# Simple figure with just one projection

# Option 1: Create a projection manually with plot.Proj()
# immport proplot as plot
# proj = plot.Proj('robin', lon_0=180)
# fig, axs = plot.subplots(nrows=2, axwidth=3, proj=proj)

# Option 2: Pass the name to 'proj' and keyword arguments to 'proj_kw'
import proplot as plot
fig, axs = plot.subplots(nrows=2, axwidth=3, proj='robin', proj_kw={'lon_0': 180})
axs.format(
    suptitle='Figure with single projection',
    coast=True, latlines=30, lonlines=60,
)

# %%
# Complex figure with different projections
import proplot as plot
fig, axs = plot.subplots(
    hratios=(1.5, 1, 1, 1, 1.5),
    basemap={
        (1, 3, 5, 7, 9): False,  # use cartopy in column 1
        (2, 4, 6, 8, 10): True,  # use basemap in column 2
    },
    proj={
        (1, 2): 'mill',  # different projection each row
        (3, 4): 'cyl',
        (5, 6): 'moll',
        (7, 8): 'sinu',
        (9, 10): 'npstere'
    },
    ncols=2, nrows=5
)
axs.format(
    suptitle='Figure with several projections',
    coast=True, latlines=30, lonlines=60,
    lonlabels='b', latlabels='r',  # or lonlabels=True, labels=True, etc.
)
axs[-1, :].format(labels=True, lonlines=30)
axs.format(collabels=['Cartopy projections', 'Basemap projections'])
plot.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geoplot:
#
# Plotting geographic data
# ------------------------
#
# In ProPlot, plotting in `~proplot.axes.GeoAxes` looks pretty much exactly the
# same as plotting in `~proplot.axes.CartesianAxes`. While cartopy and basemap
# assume your data is in "map projection" coordinates unless specified otherwise,
# ProPlot makes longitude-latitude (i.e. Plate Carrée) coordinates the *default*
# coordinate system for your datasets by passing ``transform=ccrs.PlateCarree()``
# to cartopy plotting methods and ``latlon=True`` to basemap plotting methods.
#
# There are also a couple plotting features specific to `~proplot.axes.GeoAxes`.
# To ensure a 2D plot like `~matplotlib.axes.Axes.contour` covers the entire globe,
# pass ``globe=True`` to the plotting command. This interpolates your data
# to the poles and the longitude seams before plotting.
#
# To mask out the parts of your data over the land or the ocean, toggle
# the :rcraw:`land` or and :rcraw:`ocean` settings and make sure the corresponding
# `zorder <https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
# is high enough to sit above all plotted content, e.g. with
# ``plot.rc.update({'land': True, 'land.zorder': 5})``.
# See the :ref:`next section <ug_geoformat>` for details.

# %%
import proplot as plot
import numpy as np

# Fake data with unusual longitude seam location and
# without coverage up to longitude seam and poles
offset = -40
lon = plot.arange(offset, 360 + offset - 1, 60)
lat = plot.arange(-60, 60 + 1, 30)
state = np.random.RandomState(51423)
data = state.rand(len(lat), len(lon))

# Plot data both without and with globe=True
for globe in (False, True,):
    fig, axs = plot.subplots(
        ncols=2, nrows=2, axwidth=2.5,
        proj='kav7', basemap={(1, 3): False, (2, 4): True}
    )
    for i, ax in enumerate(axs):
        cmap = ('sunset', 'sunrise')[i % 2]
        if i < 2:
            m = ax.contourf(lon, lat, data, cmap=cmap, globe=globe, extend='both')
            fig.colorbar(m, loc='b', span=i + 1, label='values', extendsize='1.7em')
        else:
            ax.pcolor(lon, lat, data, cmap=cmap, globe=globe, extend='both')
    string = 'with' if globe else 'without'
    axs.format(
        suptitle=f'Geophysical data {string} global coverage',
        collabels=['Cartopy example', 'Basemap example'],
        rowlabels=['Contourf', 'Pcolor'],
        abc=True, abcstyle='a)', abcloc='ul', abcborder=False,
        land=True, landzorder=3, lonlines=90,
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geoformat:
#
# Formatting projections
# ----------------------
#
# `~proplot.axes.CartopyAxes` and `~proplot.axes.BasemapAxes` both derive
# from `~proplot.axes.GeoAxes`, which provides the
# `proplot.axes.GeoAxes.format` method. This can
# be used to draw "major" gridlines "minor" gridlines with the `grid`,
# `longrid`, `latgrid`, and `gridminor`, `longridminor`, and `latgridminor`
# keywords. Gridline locations and label formats are configured with the
# `lonlocator`, `latlocator`, `lonformatter`, `latformatter`, `lonminorlocator`,
# and `latminorlocator` keywords. Major gridline labels and their positions
# can be configured with the `labels`, `lonlabels`, and `latlabels` keywords.
# Cartopy map bounds can be set with the `lonlim` and `latlim` keywordis. Geographic
# features like land masses, coastlines, and administrative borders can be toggled
# on and off and stylized with a variety of :ref:`rc settings <rc_proplot>`.
# Finally, `proplot.axes.GeoAxes.format` also calls `proplot.axes.Axes.format`,
# and so can be used to for subplot titles, a-b-c labels, and figure titles
# as before.
#
# For details, see the `proplot.axes.GeoAxes.format` documentation.

# %%
import proplot as plot
fig, axs = plot.subplots(
    [[1, 1, 2], [3, 3, 3]],
    axwidth=4, proj={1: 'eqearth', 2: 'ortho', 3: 'wintri'},
    wratios=(1, 1, 1.2), hratios=(1, 1.2),
)
axs.format(
    suptitle='Projection axes formatting demo',
    collabels=['Column 1', 'Column 2'],
    abc=True, abcstyle='A.', abcloc='ul', abcborder=False, linewidth=1.5
)

# Styling projections in different ways
ax = axs[0]
ax.format(
    title='Equal earth', land=True, landcolor='navy', facecolor='pale blue',
    coastcolor='gray5', borderscolor='gray5', innerborderscolor='gray5',
    gridlinewidth=1.5, gridcolor='gray5', gridalpha=0.5,
    gridminor=True, gridminorlinewidth=0.5,
    coast=True, borders=True, borderslinewidth=0.8,
)
ax = axs[1]
ax.format(
    title='Orthographic', reso='med', land=True, coast=True, latlines=10, lonlines=15,
    landcolor='mushroom', suptitle='Projection axes formatting demo',
    facecolor='petrol', coastcolor='charcoal', coastlinewidth=0.8, gridlinewidth=1
)
ax = axs[2]
ax.format(
    land=True, facecolor='ocean blue', landcolor='bisque', title='Winkel tripel',
    lonlines=60, latlines=15,
    gridlinewidth=0.8, gridminor=True, gridminorlinestyle=':',
    lonlabels=True, latlabels='r', loninline=True,
    gridlabelcolor='gray8', gridlabelsize='med-large',
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
# subintervals by passing ``dms=False`` to `~proplot.axes.GeoAxes.format`.
#
# To zoom into basemap projections, pass any of the `boundinglat`,
# `llcrnrlon`, `llcrnrlat`, `urcrnrlon`, `urcrnrlat`, `llcrnrx`, `llcrnry`,
# `urcrnrx`, `urcrnry`, `width`, or `height` keyword arguments to
# the `~proplot.constructor.Proj` constructor function either directly or via
# the `proj_kw` `~proplot.ui.subplots` keyword argument. You can also pass
# `lonlim` and `latlim` to `~proplot.constructor.Proj` and these arguments
# will be used for `llcrnrlon`, `llcrnrlat`, etc. You can not zoom into basemap
# projections with `format` after they have already been created.

# %%
import proplot as plot

# Plate Carrée map projection
plot.rc.reso = 'med'  # use higher res for zoomed in geographic features
proj = plot.Proj('cyl', lonlim=(-20, 180), latlim=(-10, 50), basemap=True)
fig, axs = plot.subplots(nrows=2, axwidth=5, proj=('cyl', proj))
axs.format(
    land=True, labels=True, lonlines=20, latlines=20,
    gridminor=True, suptitle='Zooming into projections'
)
axs[0].format(
    lonlim=(-140, 60), latlim=(-10, 50),
    labels=True, title='Cartopy example'
)
axs[1].format(title='Basemap example')

# %%
import proplot as plot

# Pole-centered map projections
proj = plot.Proj('npaeqd', boundinglat=60, basemap=True)
fig, axs = plot.subplots(ncols=2, axwidth=2.7, proj=('splaea', proj))
axs.format(
    land=True, latmax=80,  # no gridlines poleward of 80 degrees
    suptitle='Zooming into polar projections'
)
axs[0].format(boundinglat=-60, title='Cartopy example')
axs[1].format(title='Basemap example')

# %%
import proplot as plot

# Zooming in on continents
proj1 = plot.Proj('lcc', lon_0=0)  # cartopy projection
proj2 = plot.Proj('lcc', lon_0=-100, lat_0=45, width=8e6, height=8e6, basemap=True)
fig, axs = plot.subplots(ncols=2, axwidth=3, proj=(proj1, proj2))
axs.format(suptitle='Zooming into specific regions', land=True, gridminor=True)
axs[0].format(lonlim=(-20, 50), latlim=(30, 70), title='Cartopy example')
axs[1].format(lonlines=20, title='Basemap example')

# Zooming to very small scale with degree-minute-second labels
plot.rc.reso = 'hi'
fig, axs = plot.subplots(ncols=2, axwidth=2.5, proj='cyl')
axs.format(
    land=True, labels=True, gridminor=True,
    borders=True, borderscolor='white',
    suptitle='Degree-minute-second labels',
)
axs[0].format(lonlim=(-7.5, 2), latlim=(49.5, 59))
axs[1].format(lonlim=(-6, -2), latlim=(54.5, 58.5))
plot.rc.reset()


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
# Kavrisky VII projections (i.e. ``'aitoff'``, ``'hammer'``, ``'wintri'``,
# and ``'kav7'``), as well as North and South polar versions of the Azimuthal
# Equidistant, Lambert Azimuthal Equal-Area, and Gnomic projections (i.e.
# ``'npaeqd'``, ``'spaeqd'``, ``'nplaea'``, ``'splaea'``, ``'npgnom'``, and
# ``'spgnom'``), modeled after the existing `~cartopy.crs.NorthPolarStereo`
# and `~cartopy.crs.SouthPolarStereo` projections.

# %%
import proplot as plot

# Table of cartopy projections
projs = [
    'cyl', 'merc', 'mill', 'lcyl', 'tmerc',
    'robin', 'hammer', 'moll', 'kav7', 'aitoff', 'wintri', 'sinu',
    'geos', 'ortho', 'nsper', 'aea', 'eqdc', 'lcc', 'gnom',
    'npstere', 'nplaea', 'npaeqd', 'npgnom', 'igh',
    'eck1', 'eck2', 'eck3', 'eck4', 'eck5', 'eck6'
]
fig, axs = plot.subplots(ncols=3, nrows=10, width=7, proj=projs)
axs.format(
    land=True, reso='lo', labels=False,
    suptitle='Table of cartopy projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)

# %%
import proplot as plot

# Table of basemap projections
projs = [
    'cyl', 'merc', 'mill', 'cea', 'gall', 'sinu',
    'eck4', 'robin', 'moll', 'kav7', 'hammer', 'mbtfpq',
    'geos', 'ortho', 'nsper',
    'vandg', 'aea', 'eqdc', 'gnom', 'cass', 'lcc',
    'npstere', 'npaeqd', 'nplaea'
]
fig, axs = plot.subplots(ncols=3, nrows=8, basemap=True, width=7, proj=projs)
axs.format(
    land=True, labels=False,
    suptitle='Table of basemap projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)
