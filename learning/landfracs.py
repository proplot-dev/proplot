def landfracs(lonb, latb=None, lonmin=-180, 
        directory=os.environ['HOME']+'/data/', 
        file='orography.nc', varname='lsm'):
        # file='fractional_land.0.25-deg.nc', varname='data'):
        # the orography file used by ERA-Interim is actually binary land fractions (>50%), but
        # it's high enough resolution that it shouldn't matter
    """
    Retrieve land fractions on the provided graticule.
    Input options:
        1) landfracs(lonb, latb) -- include graticule vectors
        2) landfracs(lonres, latres) -- use specified latitude, longitude resolution.
        3) landfracs(res) -- use same resolution in latitude and longitude.
    """
    # file = os.environ['HOME'] + '/landfracs.nc'
    # file = os.environ['HOME'] + '/data/
    F, f = nc.ncload(directory + file, varname, 
            graticule=True, verbose=False) # want to return graticule
        # use the ERA-Interim land-sea mask at 0.125 resolution
    if latb is None: # if, for example, want 5-by-5 resolution
        latb = lonb
    # Latitudes
    try:
        iter(lonb)
    except TypeError:
        lonb = np.arange(lonmin, lonmin+360+lonb/2, lonb) # like np.arange, but includes endpoints
    else:
        if type(lonb) is not np.ndarray: lonb = np.array(lonb)
    # Longitudes
    try:
        iter(latb)
    except TypeError:
        latb = np.arange(-90, 90+latb/2, latb) # like np.arange, but includes endpoints
    else:
        if type(latb) is not np.ndarray: latb = np.array(latb)
    # And downsample
    return downsample(f.lonb, f.latb, F, lonb, latb)

