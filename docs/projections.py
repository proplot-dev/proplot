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
# .. _ug_proj:
#
# Geographic and polar plots
# ==========================
#
# ProPlot includes features for working with `polar axes
# <https://matplotlib.org/3.1.0/gallery/pie_and_polar_charts/polar_demo.html>`__
# and the `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__ and
# `basemap <https://matplotlib.org/basemap/index.html>`__ map projection
# packages. These features are optional -- installation of cartopy and basemap
# are not required.
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
# To draw polar axes, pass ``proj='polar'`` or e.g. ``proj={1:'polar'}`` to
# `~proplot.ui.subplots`. This generates a `~proplot.axes.PolarAxes`
# instance. Its `~proplot.axes.PolarAxes.format` command permits
# polar-specific modifications like changing the central radius, the zero
# azimuth location, the radial and azimuthal limits, and the positive
# azimuthal direction. This projection can also be used to make sector plots
# and annular plots.

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
    rlabelpos=180, color='gray8', ticklabelpad='1em'
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
# To turn a subplot into a plottable map projection, pass
# ``proj='name'`` or e.g. ``proj={2: 'name'}``
# (:ref:`see above <ug_proj>`) to  `~proplot.ui.subplots`
# where ``'name'`` is any valid :ref:`PROJ projection name <proj_included>`,
# or supply `proj` with a cartopy `~cartopy.crs.Projection` or basemap
# `~mpl_toolkits.basemap.Basemap` instance returned by the
# `~proplot.constructor.Proj` constructor function. Cartopy is used by
# default, but you can switch to basemap using ``basemap=True``.
# When you request map projections, `~proplot.ui.subplots` returns
# instances of `~proplot.axes.CartopyAxes` or `~proplot.axes.BasemapAxes`.
#
# * `~proplot.axes.CartopyAxes` joins the cartopy
#   `~cartopy.mpl.geoaxes.GeoAxes` class with the ProPlot
#   `~matplotlib.axes.Axes` class and adds a `~proplot.axes.GeoAxes.format`
#   command. This class includes all the normal `~cartopy.mpl.geoaxes.GeoAxes`
#   methods, and its `~proplot.axes.GeoAxes.format` method can be used to set
#   the map bounds with `~cartopy.mpl.geoaxes.GeoAxes.set_extent` and add
#   geographic features with `~cartopy.mpl.geoaxes.GeoAxes.add_feature`.
#
#   ProPlot also makes sure polar projections like `~cartopy.crs.NorthPolarStereo`
#   have circular bounds. By default, polar Gnomonic projections are bounded
#   at 30 degrees latitude, other polar projections are bounded at the equator,
#   and all other projections are given global extent using
#   `~cartopy.mpl.geoaxes.GeoAxes.set_global`. This behavior can be disabled
#   by setting :rc:`cartopy.autoextent` to ``True``, which lets cartopy automatically
#   determine the map boundaries based on plotted content.
#
# * `~proplot.axes.BasemapAxes` redirects the plot, scatter, contour,
#   contourf, pcolor, pcolormesh, quiver, streamplot, and barb methods to
#   identically named methods on the `~mpl_toolkits.basemap.Basemap` instance,
#   and provides access to `~mpl_toolkits.basemap.Basemap` geographic plotting
#   commands like `~mpl_toolkits.basemap.Basemap.fillcontinents` via
#   `~proplot.axes.GeoAxes.format`.
#
# This means with ProPlot, you no longer have to invoke verbose cartopy
# `~cartopy.crs.Projection` classes like `~cartopy.crs.LambertAzimuthalEqualArea`,
# and you never have to directly reference the `~mpl_toolkits.basemap.Basemap`
# instance -- ProPlot works with the `~mpl_toolkits.basemap.Basemap` instance
# under the hood.
#
# To make things a bit more consistent, the `~proplot.constructor.Proj` constructor
# function lets you supply native `PROJ <https://proj.org>`__ keyword names to
# the cartopy `~cartopy.crs.Projection` classes, e.g. `lon_0` instead of
# `central_longitude`. It also lets you instantiate `~mpl_toolkits.basemap.Basemap`
# projections with sensible defaults rather than raising an error when certain
# projection arguments are omitted.

# %%
import proplot as plot
plot.rc.update({
    'geogrid.linewidth': 0.5,
    'geogrid.linestyle': '-',
    'geogrid.alpha': 0.2,
})

# Simple figure with just one projection
# Option 1: Create a projection manually with plot.Proj()
# proj = plot.Proj('robin', lon_0=180)
# fig, axs = plot.subplots(ncols=2, axwidth=2.5, proj=proj)
# Option 2: Pass the name to 'proj' and keyword arguments to 'proj_kw'
fig, axs = plot.subplots(ncols=2, axwidth=2.5, proj='robin', proj_kw={'lon_0': 180})
axs.format(
    suptitle='Figure with single projection',
    coast=True, latlines=30, lonlines=60,
)

# Complex figure with different projections
fig, axs = plot.subplots(
    hratios=(1.5, 1, 1, 1, 1),
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
    coast=True, latlines=30, lonlines=60, labels=True,
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
# The below example demonstrates how to plot geographic data with ProPlot.
# It is mostly the same as cartopy, but with some new features powered by the
# `~proplot.axes.standardize_2d`, `~proplot.axes.default_transform`,
# and `~proplot.axes.default_latlon` wrappers.
#
# * For both basemap and cartopy projections, you can pass ``globe=True`` to
#   2D plotting commands to ensure global data coverage.
# * For `~proplot.axes.CartopyAxes` plotting methods,
#   ``transform=crs.PlateCarree()`` is now the default behavior. That is,
#   ProPlot assumes your data is in longitude-latitude coordinates rather than
#   map projection coordinates.
# * For `~proplot.axes.BasemapAxes` plotting methods, ``latlon=True`` is now
#   the default behavior. Again, plotting methods are now called on the *axes*
#   instead of the `~mpl_toolkits.basemap.Basemap` instance.
#
# To mask out the parts of your data over the land or the ocean, simply toggle
# the :rcraw:`land` or and :rcraw:`ocean` settings and set the corresponding
# `zorder <https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
# to a high value with e.g. ``plot.rc.update({'land': True, 'land.zorder': 5})``.
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
titles = ('Geophysical data demo', 'Global coverage demo')
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
    axs.format(
        suptitle=titles[globe],
        collabels=['Cartopy example', 'Basemap example'],
        rowlabels=['Contourf', 'Pcolor'],
        abc=True, abcstyle='a)', abcloc='ul', abcborder=False,
        coast=True, lonlines=90,
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_geoformat:
#
# Formatting projections
# ----------------------
#
# `~proplot.axes.CartopyAxes` and `~proplot.axes.BasemapAxes` both derive
# from `~proplot.axes.GeoAxes`, which provides a
# `~proplot.axes.GeoAxes.format` method. `~proplot.axes.GeoAxes.format` can
# be used to draw "major" gridlines, "minor" gridlines, add gridline labels
# with optional degree-minute-second units (cartopy > 0.18), specify gridline
# label locations, modify the projection bounding box, and add and stylize common
# geographic features like land masses, coastlines, and administrative
# borders. This method also calls `format` on `~proplot.axes.Axes`, and so
# can be used for subplot titles, a-b-c labels, and figure titles as before.

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
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_zoom:
#
# Zooming into projections
# ------------------------
#
# To zoom into cartopy projections, you can use
# `~cartopy.mpl.geoaxes.GeoAxes.set_extent`, or alternatively pass `lonlim`,
# `latlim`, or `boundinglat` to `~proplot.axes.GeoAxes.format`. The `boundinglat`
# keyword controls the circular latitude boundary for North Polar and
# South Polar Stereographic, Azimuthal Equidistant, Lambert Azimuthal
# Equal-Area, and Gnomonic projections.
#
# To zoom into basemap projections, you can pass any of the `boundinglat`,
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
proj = plot.Proj('cyl', llcrnrlon=-20, llcrnrlat=-10, urcrnrlon=180, urcrnrlat=50, basemap=True)  # noqa: E501
fig, axs = plot.subplots(nrows=2, axwidth=4.5, proj=('cyl', proj))
axs.format(
    land=True, labels=True, lonlines=20,
    latlines=20, suptitle='Zooming into projections'
)
axs[0].format(
    lonlim=(-140, 60), latlim=(-10, 50),
    labels=True, title='Cartopy example'
)
axs[1].format(title='Basemap example')

# Pole-centered map projections
proj = plot.Proj('npaeqd', boundinglat=60, basemap=True)
fig, axs = plot.subplots(ncols=2, axwidth=2.2, proj=('splaea', proj))
axs.format(
    land=True, latlines=10, latmax=80,
    suptitle='Zooming into polar projections'
)
axs[0].format(boundinglat=-60, title='Cartopy example')
axs[1].format(title='Basemap example')

# Focusing on continents
proj1 = plot.Proj('lcc', lon_0=0)  # cartopy projection
proj2 = plot.Proj('lcc', lon_0=-100, lat_0=45, width=8e6, height=8e6, basemap=True)
fig, axs = plot.subplots(ncols=2, axwidth=2.4, proj=(proj1, proj2))
axs.format(suptitle='Zooming into specific regions', land=True)
axs[0].format(
    title='Cartopy example', land=True,
    lonlim=(-20, 50), latlim=(30, 70)
)
axs[1].format(title='Basemap example', land=True)
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
import numpy as np

# Table of cartopy projections
projs = [
    'cyl', 'merc', 'mill', 'lcyl', 'tmerc',
    'robin', 'hammer', 'moll', 'kav7', 'aitoff', 'wintri', 'sinu',
    'geos', 'ortho', 'nsper', 'aea', 'eqdc', 'lcc', 'gnom',
    'npstere', 'nplaea', 'npaeqd', 'npgnom', 'igh',
    'eck1', 'eck2', 'eck3', 'eck4', 'eck5', 'eck6'
]
fig, axs = plot.subplots(ncols=3, nrows=10, proj=projs)
axs.format(
    land=True, reso='lo', labels=False,
    suptitle='Table of cartopy projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)

# Table of basemap projections
projs = [
    'cyl', 'merc', 'mill', 'cea', 'gall', 'sinu',
    'eck4', 'robin', 'moll', 'kav7', 'hammer', 'mbtfpq',
    'geos', 'ortho', 'nsper',
    'vandg', 'aea', 'eqdc', 'gnom', 'cass', 'lcc',
    'npstere', 'npaeqd', 'nplaea'
]
fig, axs = plot.subplots(ncols=3, nrows=8, basemap=True, proj=projs)
axs.format(
    land=True, labels=False,
    suptitle='Table of basemap projections'
)
for proj, ax in zip(projs, axs):
    ax.format(title=proj, titleweight='bold', labels=False)
