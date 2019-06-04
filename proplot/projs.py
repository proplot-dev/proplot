#!/usr/bin/env python3
"""
Introduces the generalized `~mpl_toolkits.basemap` `~mpl_toolkits.basemap.Basemap`
and `cartopy.crs.Projection` projection instantiator "`Proj`".
Also "registers" cartopy projections by their `PROJ.4 aliases
<https://proj4.org/operations/projections/index.html>`_ like in
`~mpl_toolkits.basemap`.

.. I'm not sure why the architects of matplotlib and related
   libraries choose to "register" certain libraries of closely-related
   classes (e.g. axis scales, projections), but not others.

Also adds pseudocylindrical `cartopy projections
<https://scitools.org.uk/cartopy/docs/latest/crs/projections.html>`_
that are currently unavailable: Hammer, Aitoff, Kavrayskiy VII,
and Winkel tripel.

Table of projections
====================
The following is a table of registered projections, their full names (with a
link to the `PROJ.4 documentation <https://proj4.org/operations/projections/index.html>`_
if it exists), and whether they are available in the cartopy and basemap packages.
Note both basemap and cartopy use `PROJ.4` as their backends.

====================================  ===========================================================================================  =========  =======
Key                                   Name                                                                                         Cartopy    Basemap
====================================  ===========================================================================================  =========  =======
``'aea'``                             `Albers Equal Area <https://proj4.org/operations/projections/aea.html>`_                     ✓          ✓
``'aeqd'``                            `Azimuthal Equidistant <https://proj4.org/operations/projections/aeqd.html>`_                ✓          ✓
``'aitoff'``                          `Aitoff <https://proj4.org/operations/projections/aitoff.html>`_                             ✓ (added)  ✗
``'cass'``                            `Cassini-Soldner <https://proj4.org/operations/projections/cass.html>`_                      ✗          ✓
``'cea'``                             `Cylindrical Equal Area <https://proj4.org/operations/projections/cea.html>`_                ✗          ✓
``'eqc'``, ```'cyl'``, ``'pcarree'``  `Cylindrical Equidistant <https://proj4.org/operations/projections/eqc.html>`_               ✓          ✓
``'eck1'``                            `Eckert I <https://proj4.org/operations/projections/eck1.html>`_                             ✓          ✗
``'eck2'``                            `Eckert II <https://proj4.org/operations/projections/eck2.html>`_                            ✓          ✗
``'eck3'``                            `Eckert III <https://proj4.org/operations/projections/eck3.html>`_                           ✓          ✗
``'eck4'``                            `Eckert IV <https://proj4.org/operations/projections/eck4.html>`_                            ✓          ✓
``'eck5'``                            `Eckert V <https://proj4.org/operations/projections/eck5.html>`_                             ✓          ✗
``'eck6'``                            `Eckert VI <https://proj4.org/operations/projections/eck6.html>`_                            ✓          ✗
``'eqdc'``                            `Equidistant Conic <https://proj4.org/operations/projections/eqdc.html>`_                    ✓          ✓
``'eqearth'``                         `Equal Earth <https://proj4.org/operations/projections/eqearth.html>`_                       ✓          ✗
``'europp'``                          Euro PP (Europe)                                                                             ✓          ✗
``'gall'``                            `Gall Stereographic Cylindrical <https://proj4.org/operations/projections/gall.html>`_       ✗          ✓
``'geos'``                            `Geostationary <https://proj4.org/operations/projections/geos.html>`_                        ✓          ✓
``'gnom'``                            `Gnomonic <https://proj4.org/operations/projections/gnom.html>`_                             ✓          ✓
``'hammer'``                          `Hammer <https://proj4.org/operations/projections/hammer.html>`_                             ✓ (added)  ✓
``'igh'``                             `Interrupted Goode Homolosine <https://proj4.org/operations/projections/igh.html>`_          ✓          ✗
``'kav7'``                            `Kavrayskiy VII <https://proj4.org/operations/projections/kav7.html>`_                       ✓ (added)  ✓
``'laea'``                            `Lambert Azimuthal Equal Area <https://proj4.org/operations/projections/laea.html>`_         ✓          ✓
``'lcc'``                             `Lambert Conformal <https://proj4.org/operations/projections/lcc.html>`_                     ✓          ✓
``'lcyl'``                            Lambert Cylindrical                                                                          ✓          ✗
``'mbtfpq'``                          `McBryde-Thomas Flat-Polar Quartic <https://proj4.org/operations/projections/mbtfpq.html>`_  ✗          ✓
``'merc'``                            `Mercator <https://proj4.org/operations/projections/merc.html>`_                             ✓          ✓
``'mill'``                            `Miller Cylindrical <https://proj4.org/operations/projections/mill.html>`_                   ✓          ✓
``'moll'``                            `Mollweide <https://proj4.org/operations/projections/moll.html>`_                            ✓          ✓
``'npaeqd'``                          North-Polar Azimuthal Equidistant                                                            ✗          ✓
``'nplaea'``                          North-Polar Lambert Azimuthal                                                                ✗          ✓
``'npstere'``                         North-Polar Stereographic                                                                    ✓          ✓
``'nsper'``                           `Near-Sided Perspective <https://proj4.org/operations/projections/nsper.html>`_              ✓          ✓
``'osni'``                            OSNI (Ireland)                                                                               ✓          ✗
``'osgb'``                            OSGB (UK)                                                                                    ✓          ✗
``'omerc'``                           `Oblique Mercator <https://proj4.org/operations/projections/omerc.html>`_                    ✗          ✓
``'ortho'``                           `Orthographic <https://proj4.org/operations/projections/ortho.html>`_                        ✓          ✓
``'poly'``                            `Polyconic <https://proj4.org/operations/projections/poly.html>`_                            ✗          ✓
``'rotpole'``                         Rotated Pole                                                                                 ✓          ✓
``'sinu'``                            `Sinusoidal <https://proj4.org/operations/projections/sinu.html>`_                           ✓          ✓
``'spaeqd'``                          South-Polar Azimuthal Equidistant                                                            ✗          ✓
``'splaea'``                          South-Polar Lambert Azimuthal                                                                ✗          ✓
``'spstere'``                         South-Polar Stereographic                                                                    ✓          ✓
``'stere'``                           `Stereographic <https://proj4.org/operations/projections/stere.html>`_                       ✓          ✓
``'tmerc'``                           `Transverse Mercator <https://proj4.org/operations/projections/tmerc.html>`_                 ✓          ✓
``'utm'``                             `Universal Transverse Mercator <https://proj4.org/operations/projections/utm.html>`_         ✓          ✗
``'vandg'``                           `van der Grinten <https://proj4.org/operations/projections/vandg.html>`_                     ✗          ✓
``'wintri'``                          `Winkel tripel <https://proj4.org/operations/projections/wintri.html>`_                      ✓ (added)  ✗
====================================  ===========================================================================================  =========  =======
"""
import numpy as np
import matplotlib.path as mpath
import warnings
from .rcmod import rc
try:
    import cartopy.crs as ccrs
    from cartopy.crs import _WarpedRectangularProjection
except ModuleNotFoundError:
    ccrs = None
    _WarpedRectangularProjection = object
# from packaging import version
# if version.parse(cartopy.__version__) < version.parse("0.13"):
#     raise RuntimeError('Require cartopy version >=0.13.') # adds set_boundary method

# Paths for cartopy projection boundaries
# WARNING: Tempting to use classmethod mpath.Path.circle, but this ends up
# failing and drawing weird polygon. Need manual approach.
# mpath.Path([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]) # rectangle
def Circle(N=100):
    """Returns a circle `~matplotlib.path.Path`. Used as the outline
    for polar stereo, aeqd, and lambert conformal projections."""
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

# Constructor function
def Proj(name, basemap=False, **kwargs):
    """
    Returns a `~mpl_toolkits.basemap.Basemap` or `cartopy.crs.Projection`
    instance.

    Parameters
    ----------
    name : str
        The projection name.
    basemap : bool, optional
        Whether to use the basemap or cartopy package. Defaults to ``False``.
    **kwargs
        Passed to the `~mpl_toolkits.basemap.Basemap` or `cartopy.crs.Projection`
        initializers.

    Returns
    -------
    proj : `~mpl_toolkits.basemap.Basemap` or `~cartopy.crs.Projection`
        The projection instance.
    aspect : float
        The aspect ratio.
    kwextra : dict
        Extra keyword args. These are still "projection" arguments but
        need to be passed to the *axes* initializer instead of the
        `~mpl_toolkits.basemap.Basemap` or `~cartopy.crs.Projection`
        initializers. So far only used by `~proplot.axes.CartopyAxes`.

    See also
    --------
    `~proplot.proj`, `CartopyAxes`, `BasemapAxes`
    """
    # Basemap
    name = name or 'cyl'
    kwextra = {}
    if basemap:
        import mpl_toolkits.basemap as mbasemap # verify package is available
        name = _basemap_cyl.get(name, name)
        kwproj = basemap_rc.get(name, {})
        kwproj.update(kwargs)
        kwproj.update({'fix_aspect': True})
        if name in _basemap_circles:
            kwproj.update({'round': True})
        reso = kwproj.pop('resolution', None) or kwproj.pop('reso', None) or 'c'
        proj = mbasemap.Basemap(projection=name, resolution=reso, **kwproj)
        aspect = (proj.urcrnrx - proj.llcrnrx) / \
                 (proj.urcrnry - proj.llcrnry)
    # Cartopy
    else:
        import cartopy.crs as ccrs # verify package is available
        kwargs = {_crs_translate.get(key, key): value for key,value in kwargs.items()}
        crs = crs_projs.get(name, None)
        if crs is None:
            raise ValueError(f'Unknown projection "{name}". Options are: {", ".join(crs_projs.keys())}.')
        for arg in ('boundinglat', 'centrallat'):
            if arg in kwargs:
                kwextra[arg] = kwargs.pop(arg)
        proj = crs(**kwargs)
        aspect = (np.diff(proj.x_limits) / \
                  np.diff(proj.y_limits))[0]
    return proj, aspect, kwextra

# Simple projections
# Inspired by source code for Mollweide implementation
class Hammer(_WarpedRectangularProjection):
    """The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__
    projection."""
    __name__ = 'hammer'
    name = 'hammer'
    """Registered projection name."""
    def __init__(self, central_longitude=0, globe=None): #, threshold=1e2):
        proj4_params = [('proj', 'hammer'), ('lon_0', central_longitude)]
        proj4_params = {'proj':'hammer', 'lon_0':central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)
    @property
    def threshold(self): # how finely to interpolate line data, etc.
        """Projection resolution."""
        return 1e4

class Aitoff(_WarpedRectangularProjection):
    """The `Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`__
    projection."""
    __name__ = 'aitoff'
    name = 'aitoff'
    """Registered projection name."""
    def __init__(self, central_longitude=0, globe=None): #, threshold=1e2):
        proj4_params = [('proj', 'aitoff'), ('lon_0', central_longitude)]
        proj4_params = {'proj':'aitoff', 'lon_0':central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)
    @property
    def threshold(self): # how finely to interpolate line data, etc.
        """Projection resolution."""
        return 1e4

class KavrayskiyVII(_WarpedRectangularProjection):
    """The `Kavrayskiy VII <https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection>`__
    projection."""
    __name__ = 'kavrayskiyVII'
    name = 'kavrayskiyVII'
    """Registered projection name."""
    def __init__(self, central_longitude=0, globe=None):
        proj4_params = [('proj', 'kav7'), ('lon_0', central_longitude)]
        super(KavrayskiyVII, self).__init__(
            proj4_params,
            central_longitude,
            globe=globe)
    @property
    def threshold(self):
        """Projection resolution."""
        return 1e4

class WinkelTripel(_WarpedRectangularProjection):
    """The `Winkel tripel (Winkel III) <https://en.wikipedia.org/wiki/Winkel_tripel_projection>`__
    projection."""
    __name__ = 'winkeltripel'
    name = 'winkeltripel'
    """Registered projection name."""
    def __init__(self, central_longitude=0, globe=None):
        proj4_params = [('proj', 'wintri'), ('lon_0', central_longitude)]
        super(WinkelTripel, self).__init__(
            proj4_params,
            central_longitude,
            globe=globe)
    @property
    def threshold(self):
        """Projection resolution."""
        return 1e4

# Basemap stuff
_basemap_circles = (
    'npstere', 'spstere', 'nplaea',
    'splaea', 'npaeqd', 'spaeqd',
    )
_basemap_cyl = { # aliases for 'cyl', that match PROJ4 name and common name
    'eqc':     'cyl',
    'pcarree': 'cyl',
    }
basemap_rc = { # note either llcrn/urcrnr args (all 4) can be specified, or width and height can be specified
    'eck4':    {'lon_0':0},
    'geos':    {'lon_0':0},
    'hammer':  {'lon_0':0},
    'moll':    {'lon_0':0},
    'kav7':    {'lon_0':0},
    'sinu':    {'lon_0':0},
    'vandg':   {'lon_0':0},
    'mbtfpq':  {'lon_0':0},
    'robin':   {'lon_0':0},
    'ortho':   {'lon_0':0, 'lat_0':0},
    'nsper':   {'lon_0':0, 'lat_0':0},
    'aea':     {'lon_0':0, 'lat_0':90, 'width':15000e3, 'height':15000e3},
    'eqdc':    {'lon_0':0, 'lat_0':90, 'width':15000e3, 'height':15000e3},
    'cass':    {'lon_0':0, 'lat_0':90, 'width':15000e3, 'height':15000e3},
    'gnom':    {'lon_0':0, 'lat_0':90, 'width':15000e3, 'height':15000e3},
    'lcc':     {'lon_0':0, 'lat_0':90, 'width':10000e3, 'height':10000e3},
    'poly':    {'lon_0':0, 'lat_0':0, 'width':10000e3, 'height':10000e3},
    'npaeqd':  {'lon_0':0, 'boundinglat':10},
    'nplaea':  {'lon_0':0, 'boundinglat':10},
    'npstere': {'lon_0':0, 'boundinglat':10},
    'spaeqd':  {'lon_0':0, 'boundinglat':-10},
    'splaea':  {'lon_0':0, 'boundinglat':-10},
    'spstere': {'lon_0':0, 'boundinglat':-10},
    'tmerc':   {'lon_0':0, 'lat_0':0, 'width':10000e3, 'height':10000e3},
    'merc':    {'llcrnrlat':-80, 'urcrnrlat':84, 'llcrnrlon':-180, 'urcrnrlon':180},
    'omerc':   {'lat_0':0, 'lon_0':0, 'lat_1':-10, 'lat_2':10, 'lon_1':0, 'lon_2':0, 'width':10000e3, 'height':10000e3},
    }
"""Default keyword args for `~mpl_toolkits.basemap.Basemap` projections.
`~mpl_toolkits.basemap` will raise an error if you don't provide them,
so ProPlot imposes some sensible default behavior."""

# Cartopy stuff
_crs_translate = { # add to this
    'lat_0':   'central_latitude',
    'lon_0':   'central_longitude',
    'lat_min': 'min_latitude',
    'lat_max': 'max_latitude',
    }
crs_projs = {}
"""Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes."""
if ccrs:
    # Custom ones, these are always present
    crs_projs = { # interpret string, create cartopy projection
      'aitoff':  Aitoff,
      'hammer':  Hammer,
      'kav7':    KavrayskiyVII,
      'wintri':  WinkelTripel,
    }
    # Builtin ones. Some of these are unavailable in older versions, so
    # we just print warning in that case.
    _unavail = []
    for _name,_class in { # interpret string, create cartopy projection
            'aea':     'AlbersEqualArea',
            'aeqd':    'AzimuthalEquidistant',
            'cyl':     'PlateCarree', # only basemap name not matching PROJ.4
            'eck1':    'EckertI',
            'eck2':    'EckertII',
            'eck3':    'EckertIII',
            'eck4':    'EckertIV',
            'eck5':    'EckertV',
            'eck6':    'EckertVI',
            'eqc':     'PlateCarree', # actual PROJ.4 name
            'eqdc':    'EquidistantConic',
            'eqearth': 'EqualEarth', # better looking Robinson; not in basemap
            'euro':    'EuroPP', # Europe; not in basemap or PROJ.4
            'geos':    'Geostationary',
            'gnom':    'Gnomonic',
            'igh':     'InterruptedGoodeHomolosine', # not in basemap
            'laea':    'LambertAzimuthalEqualArea',
            'lcc':     'LambertConformal',
            'lcyl':    'LambertCylindrical', # not in basemap or PROJ.4
            'merc':    'Mercator',
            'mill':    'Miller',
            'moll':    'Mollweide',
            'npstere': 'NorthPolarStereo', # north/south pole stuff not in PROJ.4
            'nsper':   'NearsidePerspective',
            'ortho':   'Orthographic',
            'osgb':    'OSGB', # UK; not in basemap or PROJ.4
            'osni':    'OSNI', # Ireland; not in basemap or PROJ.4
            'pcarree': 'PlateCarree', # common alt name
            'robin':   'Robinson',
            'rotpole': 'RotatedPole',
            'sinu':    'Sinusoidal',
            'spstere': 'SouthPolarStereo',
            'stere':   'Stereographic',
            'tmerc' :  'TransverseMercator',
            'utm':     'UTM', # not in basemap
            }.items():
        _class = getattr(ccrs, _class, None)
        if _class is None:
            _unavail.append(_name)
            continue
        crs_projs[_name] = _class
    if _unavail:
        warnings.warn(f'Cartopy projection(s) {", ".join(_unavail)} are unavailable. Consider updating to the latest version of cartopy.')



