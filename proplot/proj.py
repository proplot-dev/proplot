"""
Adds pseudocylindrical `cartopy projections
<https://scitools.org.uk/cartopy/docs/latest/crs/projections.html>`_
that are currently unavailable: Hammer, Aitoff, Kavrayskiy VII,
and Winkel tripel.

Also "registers" cartopy projections by their `PROJ.4 aliases
<https://proj4.org/operations/projections/index.html>`_ like in
`~mpl_toolkits.basemap`. Projections can be drawn using the `proj`
keyword arg with the `~proplot.subplots.subplots` function.

.. I'm not sure why the architects of matplotlib and related
   libraries choose to "register" certain libraries of closely-related
   classes (e.g. axis scales, projections), but not others.

The following is a table of registered projections, their full names (with a
link to the `PROJ.4 documentation <https://proj4.org/operations/projections/index.html>`_
if it exists), and whether they are available in the cartopy and basemap packages.

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
``'omerc'``                           `Oblique Mercator <https://proj4.org/operations/projections/omerc.html>`_                    ✓          ✓
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
#------------------------------------------------------------------------------#
# This module contains some custom cartopy projections, and tools for
# making stuff easier
#------------------------------------------------------------------------------#
import numpy as np
import matplotlib.path as mpath
try:
    import cartopy.crs as ccrs
    from cartopy.crs import _WarpedRectangularProjection
except ModuleNotFoundError:
    ccrs = None
    _WarpedRectangularProjection = object

# Circle path suitable for polar stereo/aeqd/lambert conformal projections
def Circle(ax, N=100):
    """
    Returns a circle `~matplotlib.path.Path`. Used as the outline
    for polar stereo, aeqd, and lambert conformal projections.
    """
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

# Constructor function
def Proj(name, **kwargs):
    """
    Returns an instance of the cartopy `~cartopy.crs.Projection` class.

    Parameters
    ----------
    name : str
        The projection name. See the table in the `~proplot.proj`
        documentation.
    **kwargs
        Passed to `~cartopy.crs.Projection`.
    """
    class_ = projs.get(name, None)
    if name is None:
        raise ValueError(f'Unknown projection "{name}". Options are: {", ".join(projs.keys())}.')
    return class_(**kwargs)

# Simple projections
# Inspired by source code for Mollweide implementation
class Hammer(_WarpedRectangularProjection):
    """
    The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__
    projection.
    """
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
    """
    The `Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`__
    projection.
    """
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
    """
    The `Kavrayskiy VII <https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection>`__
    projection.
    """
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

# TODO: Check this, but should be pretty much identical to above
class WinkelTripel(_WarpedRectangularProjection):
    """
    The `Winkel tripel (Winkel III) <https://en.wikipedia.org/wiki/Winkel_tripel_projection>`__
    projection.
    """
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

# Dictionary of names
projs = {}
"""
Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes.
"""
if ccrs:
    projs = { # interpret string, create cartopy projection
      'aea':     ccrs.AlbersEqualArea,
      'aeqd':    ccrs.AzimuthalEquidistant,
      'aitoff':  Aitoff,
      'cyl':     ccrs.PlateCarree, # only basemap name not matching PROJ.4
      'eck1':    ccrs.EckertI,
      'eck2':    ccrs.EckertII,
      'eck3':    ccrs.EckertIII,
      'eck4':    ccrs.EckertIV,
      'eck5':    ccrs.EckertV,
      'eck6':    ccrs.EckertVI,
      'eqc':     ccrs.PlateCarree, # actual PROJ.4 name
      'eqdc':    ccrs.EquidistantConic,
      'eqearth': ccrs.EqualEarth, # better looking Robinson; not in basemap
      'euro':    ccrs.EuroPP, # Europe; not in basemap or PROJ.4
      'geos':    ccrs.Geostationary,
      'gnom':    ccrs.Gnomonic,
      'hammer':  Hammer,
      'igh':     ccrs.InterruptedGoodeHomolosine, # not in basemap
      'kav7':    KavrayskiyVII,
      'laea':    ccrs.LambertAzimuthalEqualArea,
      'lcc':     ccrs.LambertConformal,
      'lcyl':    ccrs.LambertCylindrical, # not in basemap or PROJ.4
      'merc':    ccrs.Mercator,
      'mill':    ccrs.Miller,
      'moll':    ccrs.Mollweide,
      'npstere': ccrs.NorthPolarStereo, # north/south pole stuff not in PROJ.4
      'nsper':   ccrs.NearsidePerspective,
      'ortho':   ccrs.Orthographic,
      'osgb':    ccrs.OSGB, # UK; not in basemap or PROJ.4
      'osni':    ccrs.OSNI, # Ireland; not in basemap or PROJ.4
      'pcarree': ccrs.PlateCarree, # common alt name
      'robin':   ccrs.Robinson,
      'rotpole': ccrs.RotatedPole,
      'sinu':    ccrs.Sinusoidal,
      'spstere': ccrs.SouthPolarStereo,
      'stere':   ccrs.Stereographic,
      'tmerc' :  ccrs.TransverseMercator,
      'utm':     ccrs.UTM, # not in basemap
      'wintri':  WinkelTripel,
      }
