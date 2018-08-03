#------------------------------------------------------------------------------#
# This module contains some custom cartopy projections, and tools for
# making stuff easier
#------------------------------------------------------------------------------#
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath

#------------------------------------------------------------------------------#
# Simple projections
# Inspired by source code for Mollweide implementation
#------------------------------------------------------------------------------#
class Hammer(ccrs._WarpedRectangularProjection):
    __name__ = 'hammer'
    name = 'hammer'
    def __init__(self, central_longitude=0, globe=None): #, threshold=1e2):
        # self._threshold = threshold
        proj4_params = [('proj', 'hammer'), ('lon_0', central_longitude)]
        proj4_params = {'proj':'hammer', 'lon_0':central_longitude}
        super().__init__(proj4_params, central_longitude, globe=globe)
    @property
    def threshold(self): # how finely to interpolate line data, etc.
        # return self._threshold
        return 1e4

class KavrayskiyVII(ccrs._WarpedRectangularProjection):
    def __init__(self, central_longitude=0, globe=None):
        proj4_params = [('proj', 'kav7'), ('lon_0', central_longitude)]
        super(KavrayskiyVII, self).__init__(
            proj4_params,
            central_longitude,
            globe=globe)
    @property
    def threshold(self):
        return 1e4

# TODO: Check this, but should be pretty much identical to above
class WinkelTripel(ccrs._WarpedRectangularProjection):
    def __init__(self, central_longitude=0, globe=None):
        proj4_params = [('proj', 'wintri'), ('lon_0', central_longitude)]
        super(WinkelTripel, self).__init__(
            proj4_params,
            central_longitude,
            globe=globe)
    @property
    def threshold(self):
        return 1e4

#------------------------------------------------------------------------------#
# Wrappers around existing projections
#------------------------------------------------------------------------------#
# Generate 'circle' path, to be applied in axes coordinates as the new
# GeoAxes boundary.
def circle(ax, N=100):
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

# class LambertAzimuthalEqualArea(ccrs.LambertAzimuthalEqualArea)
#     def __init__(self, *args, **kwargs):
#         super(LambertAzimuthalEqualArea, self)
#
# class AzimuthalEquidistant(ccrs.AzimuthalEquidistant)
#     def __init__(self, *args, **kwargs):
#         super(AzimuthalEquidistant, self)


