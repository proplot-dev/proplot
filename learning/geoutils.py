"""
Geographic/grid related utilities.
"""
import numpy as np
from . import CONSTANTS as c
# from .basics import Dict
# __all__ = [
#     'graticule', 'downsample',
#     'seamfix', 'lonfix', 'latfix', 'geopad',
#     'landfracs', 'geomean', 'geoweights', # geographic mean, weightings
#     'haversine',
#         ]

#------------------------------------------------------------------------------
# Definitions
#------------------------------------------------------------------------------
class Properties():
    """
    For storing latitude/longitude grid properties. Assumes global grid, and
    borders halfway between each grid center.
    """
    def __init__(self, lon, lat):
        # First, guess cell widths and edges
        dlon1, dlon2, dlat1, dlat2 = lon[1]-lon[0], lon[-1]-lon[-2], lat[1]-lat[0], lat[-1]-lat[-2]
        self.latb = np.concatenate((lat[:1]-dlat1/2, (lat[1:]+lat[:-1])/2, lat[-1:]+dlat2/2))
        self.lonb = np.concatenate((lon[:1]-dlon1/2, (lon[1:]+lon[:-1])/2, lon[-1:]+dlon2/2))
        
        # Cell centers
        self.latc, self.lonc = lat.copy(), lon.copy()

        # Corrections
        # ...switch
        has90 = True if lat[-1]==90 else False
        hasm90 = True if lat[0]==-90 else False # need these switches later
        # ...use corrections for dumb grids with 'centers' at poles
        if hasm90:
            self.latc[0], self.latb[0] = -90+dlat1/4, -90
        if has90:
            self.latc[-1], self.latb[-1] = 90-dlat2/4, 90
        # ...corrected grid widths (cells half as tall near pole)
        self.dlon = self.lonb[1:]-self.lonb[:-1]
        self.dlat = self.latb[1:]-self.latb[:-1]
        if hasm90:
            self.dlat[0] /= 2
        if has90:
            self.dlat[-1] /= 2
        
        # Theta/phi coordinates
        self.phic, self.phib, self.dphi = latc*np.pi/180, latb*np.pi/180, dlat*np.pi/180
        self.thetac, self.thetab, self.dtheta = lonc*np.pi/180, lonb*np.pi/180, dlon*np.pi/180

        # Area weights (function of latitude only)
        # ...includes the latitude correction
        self.weights = dphi[None,:]*dtheta[:,None]*np.cos(phic[None,:])
        self.areas = dphi[None,:]*dtheta[:,None]*np.cos(phic[None,:])*(c.a**2)
        # areas = dphi*dtheta*np.cos(phic)*(c.a**2)
        #     # close approximation to area; cosine is extremely accurate
        # areas = areas[None,:]
        #     # make lon by lat, so broadcasting rules can apply

    def __str__(self): # what happens when calling print()
        return ("Properties of grid latitudes/longitudes.\n"
        f"Lat centers (latc): {self.latc[0]:.2f}, {self.latc[1]:.2f}, {self.latc[2]:.2f}, ... {self.latc[-1]:.2f}\n"
        f"Lat borders (latb): {self.latb[0]:.2f}, {self.latb[1]:.2f}, {self.latb[2]:.2f}, ... {self.latb[-1]:.2f}\n"
        f"Lat widths (dlat): {self.dlat[0]:.2f}, {self.dlat[1]:.2f}, {self.dlat[2]:.2f}, ... {self.dlat[-1]:.2f}\n"
        f"Lon centers (lonc): {self.lonc[0]:.2f}, {self.lonc[1]:.2f}, {self.lonc[2]:.2f}, ... {self.lonc[-1]:.2f}\n"
        f"Lon borders (lonb): {self.lonb[0]:.2f}, {self.lonb[1]:.2f}, {self.lonb[2]:.2f}, ... {self.lonb[-1]:.2f}\n"
        f"Lon centers (dlon): {self.dlon[0]:.2f}, {self.dlon[1]:.2f}, {self.dlon[2]:.2f}, ... {self.dlon[-1]:.2f}\n"
        "Coordinates in radians also provided (lat=phi, lon=theta).\n"
        "Approximate grid cell areas also provided (longitude x latitude).\n")
    def __repr__(self): # what happens when viewing in interactive session
        return self.__str__()

def pad(lon, lat, data, nlon=1, nlat=0):
    """
    Returns array padded circularly along lons, and 
    over the earth pole, for finite difference methods.
    """
    # Get padded array
    if nlon>0:
        pad = ((nlon,nlon),) + (data.ndim-1)*((0,0),)
        data = np.pad(data, pad, mode='wrap')
        lon = np.pad(lon, nlon, mode='wrap') # should be vector
    if nlat>0:
        if (data.shape[0] % 2)==1:
            raise ValueError('Data must have even number of longitudes, if you wish to pad over the poles.')
        data_append = np.roll(np.flip(data, axis=1), data.shape[0]//2, axis=0)
            # data is now descending in lat, and rolled 180 degrees in lon
        data = np.concatenate((
            data_append[:,-nlat:,...], # -87.5, -88.5, -89.5 (crossover)
            data, # -89.5, -88.5, -87.5, ..., 87.5, 88.5, 89.5 (crossover)
            data_append[:,:nlat,...], # 89.5, 88.5, 87.5
            ), axis=1)
        lat = np.pad(lat, nlat, mode='symmetric')
        lat[:nlat], lat[-nlat:] = 180-lat[:nlat], 180-lat[-nlat:]
            # much simpler for lat; but want monotonic ascent, so allow these to be >90, <-90
    return lon, lat, data

def mean(lon, lat, data, box=(None,None,None,None),
        landfracs=None, mode=None, weights=1, keepdims=False):
    """
    Takes area mean of data time series; zone and F are 2d, but data is 3d.
    Since array is masked, this is super easy... just use the masked array
    implementation of the mean, and it will adjust weights accordingly.
    
    lon: grid longitude centers
    lat: grid latitude centers
    mode: weight by land/sea coverage, or not at all
    landfracs: to be used by mode, above
    data: should be lon by lat by (extra dims); will preserve dimensionality.
    box: mean-taking region; note if edges don't fall on graticule, will just subsample
        the grid cells that do -- haven't bothered to account for partial coverage, because
        it would be pain in the butt and not that useful.
    weights: extra, custom weights to apply to data -- could be land/ocean fractions, for example.
    
    Data should be loaded with myfuncs.ncutils.ncload, and provide this function
    with metadata in 'm' structure.
    """
    # Get cell areas
    a = areas(lon, lat)

    # Get zone for average
    delta = lat[1]-lat[0]
    zone = np.ones((1,lat.size), dtype=bool)
    if box[1] is not None: # south
        zone = zone & ((lat[:-1]-delta/2 >= box[1])[None,:])
    if box[3] is not None: # north
        zone = zone & ((lat[1:]+delta/2 <= box[3])[None,:])

    # Get lune for average
    delta = lon[1]-lon[0]
    lune = np.ones((lon.size,1), dtype=bool)
    if box[0] is not None: # west
        lune = lune & ((lon[:-1]-delta/2 >= box[0])[:,None]) # left-hand box edges >= specified edge
    if box[2] is not None: # east
        lune = lune & ((lon[1:]+delta/2 <= box[2])[:,None]) # right-hand box edges <= specified edge

    # Get total weights
    weights = a*zone*lune*weights
    for s in data.shape[2:]:
        weights = weights[...,None]
    weights = np.tile(weights, (1,1,*data.shape[2:])) # has to be tiled, to match shape of data exactly

    # And apply; note that np.ma.average was tested to be slower than the 
    # three-step procedure for NaN ndarrays used below
    if type(data) is np.ma.MaskedArray:
        data = data.filled(np.nan)
    isvalid = np.isfinite(data) # boolean function extremely fast
    data[~isvalid] = 0 # scalar assignment extremely fast
    try: ave = np.average(data, weights=weights*isvalid, axis=(0,1))
    except ZeroDivisionError:  
        ave = np.nan # weights sum to zero

    # Return, optionally restoring the lon/lat dimensions to singleton
    if keepdims: 
        return ave[None,None,...]
    else:
        return ave

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    
    Input...
        lon1, lat1, lon2, lat2 (positional): ndarrays of longitude and latitude
        in degrees...
            -each **pair** should have identical shape
            -if both pairs are non-scalar, they must **also** be identically shaped
            -if one pair is scalar, distances are calculated between the scalar pair
                and each element of the non-scalar pari.
    Output...
        d: ndarray of distances in km
    """
    # Earth radius
    R = 6371.
    # Convert to radians, get differences
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    # Haversine, in km
    km = 2.*R*np.arcsin(np.sqrt(
        np.sin(dlat/2.)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.)**2
        ))
    return km

