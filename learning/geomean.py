# COPY OF GEOMEAN; USES LONB/LATB INSTEAD OF JUST CENTERS.
def mean(lon, lat, data, box=(None,None,None,None),
        landfracs=None, mode=None, weights=1, keepdims=False):
    """
    Takes area mean of data time series; zone and F are 2d, but data is 3d.
    Since array is masked, this is super easy... just use the masked array
    implementation of the mean, and it will adjust weights accordingly.
    
    lonb: longitude graticule
    latb: latitude graticule
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
    # First get lon, lat borders and box widths for weighting
    # ...use corrections for dumb grids with centers at poles
    lonb, latb = basics.graticule(lon, lat)
    if lat[0]==90:
        lat[0] = -90 + (lat[1]-lat[0])/4
        latb[0] = -90
    if lat[-1]==-90:
        lat[-1] = 90 + (lat[-1]-lat[-2])/4
        latb[-1] = 90
    # ...and next, the deltas
    delta_lon = (lonb[1:]-lonb[:-1])[...,np.newaxis]*np.pi/180
    delta_lat = (latb[1:]-latb[:-1])[np.newaxis,...]*np.pi/180
    lat = (latb[1:]+latb[:-1])/2

    # Get zone
    zone = np.ones((1,latb.size-1),dtype=bool)
    if box[1] is not None: # south
        zone = zone & ((latb[:-1]>=box[1])[np.newaxis,...])
    if box[3] is not None: # north
        zone = zone & ((latb[1:]<=box[3])[np.newaxis,...])

    # Get lune
    lune = np.ones((lonb.size-1,1),dtype=bool)
    if box[0] is not None: # west
        lune = lune & ((lonb[:-1]>=box[0])[...,np.newaxis]) # left-hand box edges >= specified edge
    if box[2] is not None: # east
        lune = lune & ((lonb[1:]<=box[2])[...,np.newaxis]) # right-hand box edges <= specified edge

    # Get area weightings (cosine is really really accurate)
    cos = np.abs(np.cos(lat[np.newaxis,...]*np.pi/180))

    # Get total weights
    w = delta_lat*delta_lon*cos*zone*lune*weights
    for s in data.shape[2:]:
        w = w[...,np.newaxis]
    w = np.tile(w, (1,1,*data.shape[2:])) # has to be tiled, to match shape of data exactly

    # And apply, allowing for different data types; will squeeze axes 0 and 1
    if type(data) is np.ma.MaskedArray:
        ave = np.ma.average(data, weights=w, axis=(0,1)) # automatically accounts for missing cells
    elif type(data) is np.ndarray: # about the same speed, actually
        filt = np.isfinite(data) # filter
        data[~filt] = 0 # this assignment is pretty fast
        try:
            ave = np.average(data, weights=w*filt, axis=(0,1)) # and the multiplication is pretty fast
        except ZeroDivisionError: # weights sum to zero
            ave = np.nan
    else:
        raise ValueError('Invalid "data" type.')

    # Optionally keep dimensions
    if keepdims: 
        ave = ave[np.newaxis,np.newaxis,...]

    # Return
    return ave
