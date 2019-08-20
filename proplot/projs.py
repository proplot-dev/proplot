#!/usr/bin/env python3
"""
Introduces the generalized `~mpl_toolkits.basemap.Basemap`
and `cartopy.crs.Projection` projection instantiator "`Proj`".
Also "registers" cartopy projections by their `PROJ.4 aliases
<https://proj4.org/operations/projections/index.html>`_ like in
`~mpl_toolkits.basemap`.

Also adds pseudocylindrical `cartopy projections
<https://scitools.org.uk/cartopy/docs/latest/crs/projections.html>`_
that are currently unavailable (Hammer, Aitoff, Kavrayskiy VII,
and Winkel tripel), along with polar versions of the Azimuthal Equidistant and
Lambert Azimuthal Equal Area projections, just like in `~mpl_toolkits.basemap`.

####################
Table of projections
####################
The following is a table of registered projections, their full names (with a
link to the `PROJ.4 documentation <https://proj4.org/operations/projections/index.html>`_
if it exists), and whether they are available in the cartopy and basemap packages.
Note both basemap and cartopy use `PROJ.4` as their backends. The ``(added)``
parenthetical indicates a cartopy projection that has been added by ProPlot.

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
``'npaeqd'``                          North-Polar Azimuthal Equidistant                                                            ✓ (added)  ✓
``'npgnom'``                          North-Polar Gnomonic                                                                         ✓ (added)  ✗
``'nplaea'``                          North-Polar Lambert Azimuthal                                                                ✓ (added)  ✓
``'npstere'``                         North-Polar Stereographic                                                                    ✓          ✓
``'nsper'``                           `Near-Sided Perspective <https://proj4.org/operations/projections/nsper.html>`_              ✓          ✓
``'osni'``                            OSNI (Ireland)                                                                               ✓          ✗
``'osgb'``                            OSGB (UK)                                                                                    ✓          ✗
``'omerc'``                           `Oblique Mercator <https://proj4.org/operations/projections/omerc.html>`_                    ✗          ✓
``'ortho'``                           `Orthographic <https://proj4.org/operations/projections/ortho.html>`_                        ✓          ✓
``'poly'``                            `Polyconic <https://proj4.org/operations/projections/poly.html>`_                            ✗          ✓
``'rotpole'``                         Rotated Pole                                                                                 ✓          ✓
``'sinu'``                            `Sinusoidal <https://proj4.org/operations/projections/sinu.html>`_                           ✓          ✓
``'spaeqd'``                          South-Polar Azimuthal Equidistant                                                            ✓ (added)  ✓
``'spgnom'``                          South-Polar Gnomonic                                                                         ✓ (added)  ✗
``'splaea'``                          South-Polar Lambert Azimuthal                                                                ✓ (added)  ✓
``'spstere'``                         South-Polar Stereographic                                                                    ✓          ✓
``'stere'``                           `Stereographic <https://proj4.org/operations/projections/stere.html>`_                       ✓          ✓
``'tmerc'``                           `Transverse Mercator <https://proj4.org/operations/projections/tmerc.html>`_                 ✓          ✓
``'utm'``                             `Universal Transverse Mercator <https://proj4.org/operations/projections/utm.html>`_         ✓          ✗
``'vandg'``                           `van der Grinten <https://proj4.org/operations/projections/vandg.html>`_                     ✗          ✓
``'wintri'``                          `Winkel tripel <https://proj4.org/operations/projections/wintri.html>`_                      ✓ (added)  ✗
====================================  ===========================================================================================  =========  =======
"""
# from packaging import version
# if version.parse(cartopy.__version__) < version.parse("0.13"):
#     raise RuntimeError('Require cartopy version >=0.13.') # adds set_boundary method
import numpy as np
import matplotlib.path as mpath
import warnings
__all__ = [
    'Circle', 'Proj',
    'basemap_rc', 'cartopy_projs',
    'Aitoff', 'Hammer', 'KavrayskiyVII',
    'NorthPolarAzimuthalEquidistant', 'NorthPolarLambertAzimuthalEqualArea',
    'SouthPolarAzimuthalEquidistant', 'SouthPolarLambertAzimuthalEqualArea',
    'WinkelTripel',
    ]
try:
    from cartopy.crs import (_WarpedRectangularProjection,
        LambertAzimuthalEqualArea, AzimuthalEquidistant, Gnomonic)
    _cartopy_installed = True
except ModuleNotFoundError:
    _WarpedRectangularProjection = object
    LambertAzimuthalEqualArea = object
    AzimuthalEquidistant = object
    Gnomonic = object
    _cartopy_installed = False

def Circle(N=100):
    """Returns a circle `~matplotlib.path.Path` used as the outline
    for polar stereographic, azimuthal equidistant, and Lambert
    conformal projections. This was developed from `this cartopy example
    <https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`_."""
    # WARNING: Tempting to use classmethod mpath.Path.circle, but this ends up
    # failing and drawing weird polygon. Need manual approach.
    # mpath.Path([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]) # rectangle
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

def Proj(name, basemap=False, **kwargs):
    """
    Returns a `~mpl_toolkits.basemap.Basemap` or `cartopy.crs.Projection`
    instance, used to interpret the `proj` and `proj_kw` arguments when
    passed to `~proplot.subplots.subplots`.

    Parameters
    ----------
    name : str
        The projection name. For a table of valid names, see the
        `~proplot.projs` documentation.
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
        The map projection aspect ratio. This may change if user zooms
        into a projection, but returning an initial guess helps avoid
        unnecessarily calculations down the line.

    See also
    --------
    `~proplot.axes.CartopyProjectionAxes`, `~proplot.axes.BasemapProjectionAxes`
    """
    # Basemap
    if basemap:
        import mpl_toolkits.basemap as mbasemap # verify package is available
        name = _basemap_cyl.get(name, name)
        kwproj = basemap_rc.get(name, {})
        kwproj.update(kwargs)
        kwproj.setdefault('fix_aspect', True)
        if name in _basemap_circles:
            kwproj.setdefault('round', True)
        if name == 'geos': # fix non-conda installed basemap issue: https://github.com/matplotlib/basemap/issues/361
            kwproj.setdefault('rsphere', (6378137.00,6356752.3142))
        reso = kwproj.pop('resolution', None) or kwproj.pop('reso', None) or 'c'
        proj = mbasemap.Basemap(projection=name, resolution=reso, **kwproj)
        aspect = (proj.urcrnrx - proj.llcrnrx) / \
                 (proj.urcrnry - proj.llcrnry)
    # Cartopy
    else:
        import cartopy.crs as ccrs # verify package is available
        kwargs = {_crs_translate.get(key, key): value for key,value in kwargs.items()}
        crs = cartopy_projs.get(name, None)
        if name == 'geos': # fix common mistake
            kwargs.pop('central_latitude', None)
        if 'boundinglat' in kwargs:
            raise ValueError(f'"boundinglat" must be passed to the ax.format() command for cartopy axes.')
        if crs is None:
            raise ValueError(f'Unknown projection "{name}". Options are: {", ".join(cartopy_projs.keys())}.')
        proj = crs(**kwargs)
        aspect = (np.diff(proj.x_limits) / \
                  np.diff(proj.y_limits))[0]
    return proj, aspect

# Various pseudo-rectangular projections
# Inspired by source code for Mollweide implementation
class Hammer(_WarpedRectangularProjection):
    """The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__
    projection."""
    __name__ = 'hammer'
    name = 'hammer'
    """Registered projection name."""
    def __init__(self, central_longitude=0, globe=None): #, threshold=1e2):
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
        proj4_params = {'proj':'kav7', 'lon_0':central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)
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
        proj4_params = {'proj':'wintri', 'lon_0':central_longitude}
        super(WinkelTripel, self).__init__(proj4_params, central_longitude, globe=globe)
    @property
    def threshold(self):
        """Projection resolution."""
        return 1e4

# Extra polar projections matching basemap's options
class NorthPolarAzimuthalEquidistant(AzimuthalEquidistant):
    """Analogous to `~cartopy.crs.NorthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=90,
                central_longitude=central_longitude, globe=globe)

class SouthPolarAzimuthalEquidistant(AzimuthalEquidistant):
    """Analogous to `~cartopy.crs.SouthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=-90,
                central_longitude=central_longitude, globe=globe)

class NorthPolarLambertAzimuthalEqualArea(LambertAzimuthalEqualArea):
    """Analogous to `~cartopy.crs.NorthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=90,
                central_longitude=central_longitude, globe=globe)

class SouthPolarLambertAzimuthalEqualArea(LambertAzimuthalEqualArea):
    """Analogous to `~cartopy.crs.SouthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=-90,
                central_longitude=central_longitude, globe=globe)

class NorthPolarGnomonic(Gnomonic):
    """Analogous to `~cartopy.crs.SouthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=90,
                central_longitude=central_longitude, globe=globe)

class SouthPolarGnomonic(Gnomonic):
    """Analogous to `~cartopy.crs.SouthPolarStereo`."""
    def __init__(self, central_longitude=0.0, globe=None):
        super().__init__(central_latitude=-90,
                central_longitude=central_longitude, globe=globe)

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
cartopy_projs = {}
"""Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes."""
if _cartopy_installed:
    # Custom ones, these are always present
    import cartopy.crs as ccrs # verify package is available
    cartopy_projs = { # interpret string, create cartopy projection
      'aitoff': Aitoff,
      'hammer': Hammer,
      'kav7':   KavrayskiyVII,
      'wintri': WinkelTripel,
      'npgnom': NorthPolarGnomonic,
      'spgnom': SouthPolarGnomonic,
      'npaeqd': NorthPolarAzimuthalEquidistant,
      'spaeqd': SouthPolarAzimuthalEquidistant,
      'nplaea': NorthPolarLambertAzimuthalEqualArea,
      'splaea': SouthPolarLambertAzimuthalEqualArea,
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
        cartopy_projs[_name] = _class
    if _unavail:
        warnings.warn(f'Cartopy projection(s) {", ".join(_unavail)} are unavailable. Consider updating to cartopy >= 0.17.0.')

