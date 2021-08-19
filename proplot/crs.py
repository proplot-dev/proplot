#!/usr/bin/env python3
"""
Additional cartopy projection classes.
"""
import warnings

from .internals import ic  # noqa: F401
from .internals import _snippet_manager

try:
    from cartopy.crs import _WarpedRectangularProjection
    from cartopy.crs import AzimuthalEquidistant, Gnomonic, LambertAzimuthalEqualArea
    from cartopy.crs import NorthPolarStereo, SouthPolarStereo  # noqa: F401 (geo.py)
except ModuleNotFoundError:
    _WarpedRectangularProjection = object
    AzimuthalEquidistant = object
    Gnomonic = object
    LambertAzimuthalEqualArea = object
    NorthPolarStereo = object
    SouthPolarStereo = object

__all__ = [
    'Aitoff',
    'Hammer',
    'KavrayskiyVII',
    'WinkelTripel',
    'NorthPolarAzimuthalEquidistant',
    'SouthPolarAzimuthalEquidistant',
    'NorthPolarGnomonic',
    'SouthPolarGnomonic',
    'NorthPolarLambertAzimuthalEqualArea',
    'SouthPolarLambertAzimuthalEqualArea',
]

_reso_docstring = """
The projection resolution.
"""
_init_docstring = """
Parameters
----------
central_longitude : float, optional
    The central meridian longitude in degrees. Default is 0.
false_easting: float, optional
    X offset from planar origin in metres. Default is 0.
false_northing: float, optional
    Y offset from planar origin in metres. Default is 0.
globe : `~cartopy.crs.Globe`, optional
    If omitted, a default globe is created.
"""
_snippet_manager['proj.reso'] = _reso_docstring
_snippet_manager['proj.init'] = _init_docstring


class Aitoff(_WarpedRectangularProjection):
    """
    The `Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`__ projection.
    """
    #: Registered projection name.
    name = 'aitoff'

    @_snippet_manager
    def __init__(
        self, central_longitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        %(proj.init)s
        """
        from cartopy._crs import Globe
        from cartopy.crs import WGS84_SEMIMAJOR_AXIS
        if globe is None:
            globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS, ellipse=None)

        a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
        b = globe.semiminor_axis or a
        if b != a or globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle elliptical globes.'
            )

        proj4_params = {'proj': 'aitoff', 'lon_0': central_longitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @_snippet_manager
    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """
        %(proj.reso)s
        """
        return 1e5


class Hammer(_WarpedRectangularProjection):
    """
    The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__ projection.
    """
    #: Registered projection name.
    name = 'hammer'

    @_snippet_manager
    def __init__(
        self, central_longitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        %(proj.init)s
        """
        from cartopy._crs import Globe
        from cartopy.crs import WGS84_SEMIMAJOR_AXIS
        if globe is None:
            globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS, ellipse=None)

        a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
        b = globe.semiminor_axis or a
        if b != a or globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle elliptical globes.'
            )

        proj4_params = {'proj': 'hammer', 'lon_0': central_longitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @_snippet_manager
    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """
        %(proj.reso)s
        """
        return 1e5


class KavrayskiyVII(_WarpedRectangularProjection):
    """
    The `Kavrayskiy VII \
<https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection>`__ projection.
    """
    #: Registered projection name.
    name = 'kavrayskiyVII'

    @_snippet_manager
    def __init__(
        self, central_longitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        %(proj.init)s
        """
        from cartopy._crs import Globe
        from cartopy.crs import WGS84_SEMIMAJOR_AXIS
        if globe is None:
            globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS, ellipse=None)

        a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
        b = globe.semiminor_axis or a
        if b != a or globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle elliptical globes.'
            )

        proj4_params = {'proj': 'kav7', 'lon_0': central_longitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @_snippet_manager
    @property
    def threshold(self):
        """
        %(proj.reso)s
        """
        return 1e5


class WinkelTripel(_WarpedRectangularProjection):
    """
    The `Winkel tripel (Winkel III) \
<https://en.wikipedia.org/wiki/Winkel_tripel_projection>`__ projection.
    """
    #: Registered projection name.
    name = 'winkeltripel'

    @_snippet_manager
    def __init__(
        self, central_longitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        %(proj.init)s
        """
        from cartopy._crs import Globe
        from cartopy.crs import WGS84_SEMIMAJOR_AXIS
        if globe is None:
            globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS, ellipse=None)

        a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
        b = globe.semiminor_axis or a
        if b != a or globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle '
                'elliptical globes.'
            )

        proj4_params = {'proj': 'wintri', 'lon_0': central_longitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @_snippet_manager
    @property
    def threshold(self):
        """
        %(proj.reso)s
        """
        return 1e5


class NorthPolarAzimuthalEquidistant(AzimuthalEquidistant):
    """
    Analogous to `~cartopy.crs.NorthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=90,
            central_longitude=central_longitude, globe=globe
        )


class SouthPolarAzimuthalEquidistant(AzimuthalEquidistant):
    """
    Analogous to `~cartopy.crs.SouthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=-90,
            central_longitude=central_longitude, globe=globe
        )


class NorthPolarLambertAzimuthalEqualArea(LambertAzimuthalEqualArea):
    """
    Analogous to `~cartopy.crs.NorthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=90,
            central_longitude=central_longitude, globe=globe
        )


class SouthPolarLambertAzimuthalEqualArea(LambertAzimuthalEqualArea):
    """
    Analogous to `~cartopy.crs.SouthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=-90,
            central_longitude=central_longitude, globe=globe
        )


class NorthPolarGnomonic(Gnomonic):
    """
    Analogous to `~cartopy.crs.NorthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=90,
            central_longitude=central_longitude, globe=globe
        )


class SouthPolarGnomonic(Gnomonic):
    """
    Analogous to `~cartopy.crs.SouthPolarStereo`.
    """
    @_snippet_manager
    def __init__(self, central_longitude=0.0, globe=None):
        """
        %(proj.init)s
        """
        super().__init__(
            central_latitude=-90,
            central_longitude=central_longitude, globe=globe
        )
