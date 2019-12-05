#!/usr/bin/env python3
"""
New cartopy projection classes and related tools.
Includes projection constructor function for generating
`~mpl_toolkits.basemap.Basemap` and cartopy `~cartopy.crs.Projection` classes
with their `PROJ.4 <https://proj4.org/operations/projections/index.html>`__
string name aliases, just like `~mpl_toolkits.basemap`.
"""
# from packaging import version
# if version.parse(cartopy.__version__) < version.parse("0.13"):
# raise RuntimeError('Require cartopy version >=0.13.') # adds
# set_boundary method
import numpy as np
import warnings
__all__ = [
    'Proj',
    'basemap_rc', 'cartopy_projs',
    'Aitoff', 'Hammer', 'KavrayskiyVII',
    'NorthPolarAzimuthalEquidistant', 'NorthPolarLambertAzimuthalEqualArea',
    'SouthPolarAzimuthalEquidistant', 'SouthPolarLambertAzimuthalEqualArea',
    'WinkelTripel',
]
try:
    from cartopy.crs import (
        _WarpedRectangularProjection,
        LambertAzimuthalEqualArea, AzimuthalEquidistant, Gnomonic)
    _cartopy_installed = True
except ModuleNotFoundError:
    _WarpedRectangularProjection = object
    LambertAzimuthalEqualArea = object
    AzimuthalEquidistant = object
    Gnomonic = object
    _cartopy_installed = False


def Proj(name, basemap=False, **kwargs):
    """
    Returns a `~mpl_toolkits.basemap.Basemap` or `cartopy.crs.Projection`
    instance, used to interpret the `proj` and `proj_kw` arguments when
    passed to `~proplot.subplots.subplots`.

    Parameters
    ----------
    name : str
        The projection name. Like basemap, we use the PROJ.4 shorthands.

        The following table lists the valid projection names, their full names
        (with links to the relevant `PROJ.4 documentation \
<https://proj4.org/operations/projections/index.html>`__),
        and whether they are available in the cartopy and basemap packages.

        ``(added)`` indicates a projection class that ProPlot has "added"
        to cartopy using the cartopy API.

        =============  ============================================================================================  =========  =======
        Key            Name                                                                                          Cartopy    Basemap
        =============  ============================================================================================  =========  =======
        ``'aea'``      `Albers Equal Area <https://proj4.org/operations/projections/aea.html>`__                     ✓          ✓
        ``'aeqd'``     `Azimuthal Equidistant <https://proj4.org/operations/projections/aeqd.html>`__                ✓          ✓
        ``'aitoff'``   `Aitoff <https://proj4.org/operations/projections/aitoff.html>`__                             ✓ (added)  ✗
        ``'cass'``     `Cassini-Soldner <https://proj4.org/operations/projections/cass.html>`__                      ✗          ✓
        ``'cea'``      `Cylindrical Equal Area <https://proj4.org/operations/projections/cea.html>`__                ✗          ✓
        ``'cyl'``      `Cylindrical Equidistant <https://proj4.org/operations/projections/eqc.html>`__               ✓          ✓
        ``'eck1'``     `Eckert I <https://proj4.org/operations/projections/eck1.html>`__                             ✓          ✗
        ``'eck2'``     `Eckert II <https://proj4.org/operations/projections/eck2.html>`__                            ✓          ✗
        ``'eck3'``     `Eckert III <https://proj4.org/operations/projections/eck3.html>`__                           ✓          ✗
        ``'eck4'``     `Eckert IV <https://proj4.org/operations/projections/eck4.html>`__                            ✓          ✓
        ``'eck5'``     `Eckert V <https://proj4.org/operations/projections/eck5.html>`__                             ✓          ✗
        ``'eck6'``     `Eckert VI <https://proj4.org/operations/projections/eck6.html>`__                            ✓          ✗
        ``'eqdc'``     `Equidistant Conic <https://proj4.org/operations/projections/eqdc.html>`__                    ✓          ✓
        ``'eqc'``      `Cylindrical Equidistant <https://proj4.org/operations/projections/eqc.html>`__               ✓          ✓
        ``'eqearth'``  `Equal Earth <https://proj4.org/operations/projections/eqearth.html>`__                       ✓          ✗
        ``'europp'``   Euro PP (Europe)                                                                              ✓          ✗
        ``'gall'``     `Gall Stereographic Cylindrical <https://proj4.org/operations/projections/gall.html>`__       ✗          ✓
        ``'geos'``     `Geostationary <https://proj4.org/operations/projections/geos.html>`__                        ✓          ✓
        ``'gnom'``     `Gnomonic <https://proj4.org/operations/projections/gnom.html>`__                             ✓          ✓
        ``'hammer'``   `Hammer <https://proj4.org/operations/projections/hammer.html>`__                             ✓ (added)  ✓
        ``'igh'``      `Interrupted Goode Homolosine <https://proj4.org/operations/projections/igh.html>`__          ✓          ✗
        ``'kav7'``     `Kavrayskiy VII <https://proj4.org/operations/projections/kav7.html>`__                       ✓ (added)  ✓
        ``'laea'``     `Lambert Azimuthal Equal Area <https://proj4.org/operations/projections/laea.html>`__         ✓          ✓
        ``'lcc'``      `Lambert Conformal <https://proj4.org/operations/projections/lcc.html>`__                     ✓          ✓
        ``'lcyl'``     Lambert Cylindrical                                                                           ✓          ✗
        ``'mbtfpq'``   `McBryde-Thomas Flat-Polar Quartic <https://proj4.org/operations/projections/mbtfpq.html>`__  ✗          ✓
        ``'merc'``     `Mercator <https://proj4.org/operations/projections/merc.html>`__                             ✓          ✓
        ``'mill'``     `Miller Cylindrical <https://proj4.org/operations/projections/mill.html>`__                   ✓          ✓
        ``'moll'``     `Mollweide <https://proj4.org/operations/projections/moll.html>`__                            ✓          ✓
        ``'npaeqd'``   North-Polar Azimuthal Equidistant                                                             ✓ (added)  ✓
        ``'npgnom'``   North-Polar Gnomonic                                                                          ✓ (added)  ✗
        ``'nplaea'``   North-Polar Lambert Azimuthal                                                                 ✓ (added)  ✓
        ``'npstere'``  North-Polar Stereographic                                                                     ✓          ✓
        ``'nsper'``    `Near-Sided Perspective <https://proj4.org/operations/projections/nsper.html>`__              ✓          ✓
        ``'osni'``     OSNI (Ireland)                                                                                ✓          ✗
        ``'osgb'``     OSGB (UK)                                                                                     ✓          ✗
        ``'omerc'``    `Oblique Mercator <https://proj4.org/operations/projections/omerc.html>`__                    ✗          ✓
        ``'ortho'``    `Orthographic <https://proj4.org/operations/projections/ortho.html>`__                        ✓          ✓
        ``'pcarree'``  `Cylindrical Equidistant <https://proj4.org/operations/projections/eqc.html>`__               ✓          ✓
        ``'poly'``     `Polyconic <https://proj4.org/operations/projections/poly.html>`__                            ✗          ✓
        ``'rotpole'``  Rotated Pole                                                                                  ✓          ✓
        ``'sinu'``     `Sinusoidal <https://proj4.org/operations/projections/sinu.html>`__                           ✓          ✓
        ``'spaeqd'``   South-Polar Azimuthal Equidistant                                                             ✓ (added)  ✓
        ``'spgnom'``   South-Polar Gnomonic                                                                          ✓ (added)  ✗
        ``'splaea'``   South-Polar Lambert Azimuthal                                                                 ✓ (added)  ✓
        ``'spstere'``  South-Polar Stereographic                                                                     ✓          ✓
        ``'stere'``    `Stereographic <https://proj4.org/operations/projections/stere.html>`__                       ✓          ✓
        ``'tmerc'``    `Transverse Mercator <https://proj4.org/operations/projections/tmerc.html>`__                 ✓          ✓
        ``'utm'``      `Universal Transverse Mercator <https://proj4.org/operations/projections/utm.html>`__         ✓          ✗
        ``'vandg'``    `van der Grinten <https://proj4.org/operations/projections/vandg.html>`__                     ✗          ✓
        ``'wintri'``   `Winkel tripel <https://proj4.org/operations/projections/wintri.html>`__                      ✓ (added)  ✗
        =============  ============================================================================================  =========  =======

    basemap : bool, optional
        Whether to use the basemap or cartopy package. Default is
        ``False``.
    **kwargs
        Passed to the `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class.

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
    `~proplot.axes.GeoAxes`, `~proplot.axes.BasemapAxes`
    """  # noqa
    # Basemap
    if basemap:
        import mpl_toolkits.basemap as mbasemap
        name = BASEMAP_TRANSLATE.get(name, name)
        kwproj = basemap_rc.get(name, {})
        kwproj.update(kwargs)
        kwproj.setdefault('fix_aspect', True)
        if name[:2] in ('np', 'sp'):
            kwproj.setdefault('round', True)
        # Fix non-conda installed basemap issue:
        # https://github.com/matplotlib/basemap/issues/361
        if name == 'geos':
            kwproj.setdefault('rsphere', (6378137.00, 6356752.3142))
        reso = kwproj.pop('resolution', None) or kwproj.pop(
            'reso', None) or 'c'
        proj = mbasemap.Basemap(projection=name, resolution=reso, **kwproj)
        aspect = (proj.urcrnrx - proj.llcrnrx) / \
                 (proj.urcrnry - proj.llcrnry)
    # Cartopy
    else:
        import cartopy.crs as _  # noqa
        kwargs = {CARTOPY_CRS_TRANSLATE.get(
            key, key): value for key, value in kwargs.items()}
        crs = cartopy_projs.get(name, None)
        if name == 'geos':  # fix common mistake
            kwargs.pop('central_latitude', None)
        if 'boundinglat' in kwargs:
            raise ValueError(
                f'"boundinglat" must be passed to the ax.format() command '
                'for cartopy axes.')
        if crs is None:
            raise ValueError(
                f'Unknown projection {name!r}. Options are: '
                + ', '.join(map(repr, cartopy_projs.keys())))
        proj = crs(**kwargs)
        aspect = (np.diff(proj.x_limits) / np.diff(proj.y_limits))[0]
    return proj, aspect


class Hammer(_WarpedRectangularProjection):
    """The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__
    projection."""
    __name__ = 'hammer'
    #: Registered projection name.
    name = 'hammer'

    def __init__(self, central_longitude=0, globe=None):  # , threshold=1e2):
        proj4_params = {'proj': 'hammer', 'lon_0': central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)

    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """Projection resolution."""
        return 1e4


class Aitoff(_WarpedRectangularProjection):
    """The `Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`__
    projection."""
    __name__ = 'aitoff'
    #: Registered projection name.
    name = 'aitoff'

    def __init__(self, central_longitude=0, globe=None):  # , threshold=1e2):
        proj4_params = {'proj': 'aitoff', 'lon_0': central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)

    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """Projection resolution."""
        return 1e4


class KavrayskiyVII(_WarpedRectangularProjection):
    """The `Kavrayskiy VII
    <https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection>`__ projection."""
    __name__ = 'kavrayskiyVII'
    #: Registered projection name.
    name = 'kavrayskiyVII'

    def __init__(self, central_longitude=0, globe=None):
        proj4_params = {'proj': 'kav7', 'lon_0': central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)

    @property
    def threshold(self):
        """Projection resolution."""
        return 1e4


class WinkelTripel(_WarpedRectangularProjection):
    """The `Winkel tripel (Winkel III)
    <https://en.wikipedia.org/wiki/Winkel_tripel_projection>`__ projection."""
    __name__ = 'winkeltripel'
    #: Registered projection name.
    name = 'winkeltripel'

    def __init__(self, central_longitude=0, globe=None):
        proj4_params = {'proj': 'wintri', 'lon_0': central_longitude}
        super(WinkelTripel, self).__init__(
            proj4_params, central_longitude, globe=globe)

    @property
    def threshold(self):
        """Projection resolution."""
        return 1e4


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


# Hidden constants
BASEMAP_TRANSLATE = {
    'eqc': 'cyl',
    'pcarree': 'cyl',
}
CARTOPY_CRS_TRANSLATE = {  # add to this
    'lat_0': 'central_latitude',
    'lon_0': 'central_longitude',
    'lat_min': 'min_latitude',
    'lat_max': 'max_latitude',
}

#: Default keyword args for `~mpl_toolkits.basemap.Basemap` projections.
#: `~mpl_toolkits.basemap` will raise an error if you don't provide them,
#: so ProPlot imposes some sensible default behavior.
basemap_rc = {
    'eck4': {'lon_0': 0},
    'geos': {'lon_0': 0},
    'hammer': {'lon_0': 0},
    'moll': {'lon_0': 0},
    'kav7': {'lon_0': 0},
    'sinu': {'lon_0': 0},
    'vandg': {'lon_0': 0},
    'mbtfpq': {'lon_0': 0},
    'robin': {'lon_0': 0},
    'ortho': {'lon_0': 0, 'lat_0': 0},
    'nsper': {'lon_0': 0, 'lat_0': 0},
    'aea': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'eqdc': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'cass': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'gnom': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'lcc': {'lon_0': 0, 'lat_0': 90, 'width': 10000e3, 'height': 10000e3},
    'poly': {'lon_0': 0, 'lat_0': 0, 'width': 10000e3, 'height': 10000e3},
    'npaeqd': {'lon_0': 0, 'boundinglat': 10},
    'nplaea': {'lon_0': 0, 'boundinglat': 10},
    'npstere': {'lon_0': 0, 'boundinglat': 10},
    'spaeqd': {'lon_0': 0, 'boundinglat': -10},
    'splaea': {'lon_0': 0, 'boundinglat': -10},
    'spstere': {'lon_0': 0, 'boundinglat': -10},
    'tmerc': {'lon_0': 0, 'lat_0': 0, 'width': 10000e3, 'height': 10000e3},
    'merc': {'llcrnrlat': -80, 'urcrnrlat': 84,
             'llcrnrlon': -180, 'urcrnrlon': 180},
    'omerc': {'lat_0': 0, 'lon_0': 0, 'lat_1': -10, 'lat_2': 10,
              'lon_1': 0, 'lon_2': 0, 'width': 10000e3, 'height': 10000e3},
}

#: Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes.
cartopy_projs = {}
if _cartopy_installed:
    # Custom ones, these are always present
    import cartopy.crs as ccrs  # verify package is available
    cartopy_projs = {  # interpret string, create cartopy projection
        'aitoff': Aitoff,
        'hammer': Hammer,
        'kav7': KavrayskiyVII,
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
    for _name, _class in {  # interpret string, create cartopy projection
            'aea': 'AlbersEqualArea',
            'aeqd': 'AzimuthalEquidistant',
            'cyl': 'PlateCarree',  # only basemap name not matching PROJ.4
            'eck1': 'EckertI',
            'eck2': 'EckertII',
            'eck3': 'EckertIII',
            'eck4': 'EckertIV',
            'eck5': 'EckertV',
            'eck6': 'EckertVI',
            'eqc': 'PlateCarree',  # actual PROJ.4 name
            'eqdc': 'EquidistantConic',
            'eqearth': 'EqualEarth',  # better looking Robinson; not in basemap
            'euro': 'EuroPP',  # Europe; not in basemap or PROJ.4
            'geos': 'Geostationary',
            'gnom': 'Gnomonic',
            'igh': 'InterruptedGoodeHomolosine',  # not in basemap
            'laea': 'LambertAzimuthalEqualArea',
            'lcc': 'LambertConformal',
            'lcyl': 'LambertCylindrical',  # not in basemap or PROJ.4
            'merc': 'Mercator',
            'mill': 'Miller',
            'moll': 'Mollweide',
            'npstere': 'NorthPolarStereo',  # np/sp stuff not in PROJ.4
            'nsper': 'NearsidePerspective',
            'ortho': 'Orthographic',
            'osgb': 'OSGB',  # UK; not in basemap or PROJ.4
            'osni': 'OSNI',  # Ireland; not in basemap or PROJ.4
            'pcarree': 'PlateCarree',  # common alternate name
            'robin': 'Robinson',
            'rotpole': 'RotatedPole',
            'sinu': 'Sinusoidal',
            'spstere': 'SouthPolarStereo',
            'stere': 'Stereographic',
            'tmerc': 'TransverseMercator',
            'utm': 'UTM',  # not in basemap
    }.items():
        _class = getattr(ccrs, _class, None)
        if _class is None:
            _unavail.append(_name)
            continue
        cartopy_projs[_name] = _class
    if _unavail:
        warnings.warn(
            f'Cartopy projection(s) {", ".join(map(repr, _unavail))} are '
            f'unavailable. Consider updating to cartopy >= 0.17.0.')
