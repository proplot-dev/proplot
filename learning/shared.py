# Import statement
import numpy as np

def lonfix(lon, data=None, lon_0=-180, roll=True, return_roll=False):
    """
    Shifts longitudes to the breakpoint specified by lon_0.
    Copied from block in ncutils.
    """
    if data is not None:
        if data.shape[0]!=lon.size:
            raise ValueError('Data must have longitude on first dimension.')
    # Roll longitudes so that leftmost lon center is as close as possible
    # to input lon_0, mod 360
    lonmod = lon % 360
    lon_0mod = lon_0 % 360
    lon_roll = -np.where((lonmod - lon_0mod) % 360 == 
            ((lonmod - lon_0mod) % 360).min())[0][0]
        # ...are some edge cases where might have more than one matching location; e.g.
        # if user supplied *borders*
    lon = np.roll(lon,lon_roll)
    if data is not None:
        data = np.roll(data,lon_roll,axis=0)
    # Enforce lons monotonically increasing from lon_0
    lon_changeover = np.where(np.diff(lon)<0)[0]
    if lon_changeover.size!=0:
        lon[lon_changeover[0]+1:] = lon[lon_changeover[0]+1:]+360
    lon -= np.array([360*round((lon[0]-lon_0)/360)],dtype=lon.dtype) # make lon_0 match lon[0]
    # Roll back, if roll is False; just wanted to adjust sizes; needed to roll before
    # because numpy can't do circular indexing
    if not roll:
        lon = np.roll(lon,-lon_roll)
        if data is not None:
            data = np.roll(data,-lon_roll,axis=0)
        lon_roll = 0
    # Return stuff, and remember rules of tuple expansion
    if data is not None and return_roll:
        return lon, data, lon_roll # everything
    elif data is not None:
        return lon, data # just new longitudes/data
    elif return_roll:
        return lon, lon_roll # longitudes/longitude rolling
    else:
        return lon # just new longitudes

def downsample(lonb1, latb1, data1, lonb2, latb2):
    """
    Downsamples data onto a more sparse grid using area-weighted averaging. User
    should provide the data longitude/latitude cell BOUNDARIES.
    Input...
        lonb1: cell graticule longitudes
        latb1: cell graticule latitudes
        data: data lying on above graticule, shaped lon by lat
        lonb2: interpolant; new cell graticule longitudes
        latb2: interpolant; new cell graticule latitudes
    Breakpoint of longitudes DOES NOT MATTER, as long as the graticule is circular
    i.e. the first matches the last mod 360, and longitudes are monotonic. Everything
    will be np.roll'ed to match.
    Latitude graticules must also be monotonic.
    """
     
    #--------------------------------------------------------------------------
    # Handle input
    #--------------------------------------------------------------------------
    # Longitudes, latitudes enforce float
    lonb1, latb1, lonb2, latb2 = lonb1.astype(float), latb1.astype(float), lonb2.astype(float), latb2.astype(float)
    # Fill masked array; for downsampling operation, we DO NOT want to ignore missing values,
    # so filling with np.nan first makes more sense (otherwise, operations on masked arrays are 
    # actually slightly faster than np.nansum, np.nanmean, etc.)
    if type(data1) is np.ma.MaskedArray:
        flag_masked = True
        data1 = data1.filled(np.nan)
    else:
        flag_masked = False
    # Raise errors if graticules don't cover globe
    for l in (latb1, latb2):
        if l[0]!=-90. or l[-1]!=90.:
            raise ValueError('Both latitude graticules must start at -90 and '
                    + 'end at 90 latitude.')
    for l in (lonb1, lonb2):
        if (l[0]%360) != (l[-1]%360):
            raise ValueError('Both longitude graticules must be circular.')
     
    #--------------------------------------------------------------------------
    # Fix longitude indexing for input/output
    #--------------------------------------------------------------------------
    # Since lonfix can't deal with repeating values, pretend left-hand borders
    # are grid centers. Then, add in the ending longitude.
    lonb1, data1 = lonfix(lonb1[:-1], data1, lon_0=-180.)
    lonb1 = np.concatenate((lonb1, lonb1[:1]+360.))
    lonb2, lon_roll = lonfix(lonb2[:-1], lon_0=-180., return_roll=True)
    lonb2 = np.concatenate((lonb2, lonb2[:1]+360.))
    if data1.shape[0] != lonb1.size-1:
        raise ValueError('Longitude must be on first dimension of data.')
    if data1.shape[1] != latb1.size-1:
        raise ValueError('Latitude must be on second dimension of data.')

    #--------------------------------------------------------------------------
    # Check input lonb1/latb1 graticules
    #--------------------------------------------------------------------------
    # Enforce each to have boundary exactly on -180, 180; augment
    # data grid if necessary (will raise flag to put data back later)
    # ...first, the starting lons
    if lonb1[0] != -180.:
        lonb1 = np.concatenate((np.array([-180.]),lonb1))
        lonb1[-1] = 180. # we've added a new tiny cell to the end
        pad_tuple = ((1,0),) + ((0,0),)*(data1.ndim-1)
        data1 = np.pad(data1, pad_tuple, 'wrap')
    # ...and then, the ending lons
    if lonb2[0] != -180.:
        lonb2 = np.concatenate((np.array([-180.]),lonb2))
        lonb2[-1] = 180.
        flag_cutoff = True
    else:
        flag_cutoff = False
    
    #--------------------------------------------------------------------------
    # Combo graticule
    #--------------------------------------------------------------------------
    # Build up a "combo" graticule; composed of new and old gridlines on top
    # of each other
    lon, lat = (
            np.unique(np.sort(np.concatenate((lonb1, lonb2)))), 
            np.unique(np.sort(np.concatenate((latb1, latb2))))
            )
    # Iterate through, and label which cells in combo graticule correspond to
    # which cells in old/new graticules
    lon1_idx, lat1_idx, lon2_idx, lat2_idx = (
            np.zeros(lon.size,dtype=int), np.zeros(lat.size,dtype=int),
            np.zeros(lon.size,dtype=int), np.zeros(lat.size,dtype=int)
            ) # ...and note bottom/left edges should correspond to index 0
    for i,l in enumerate(lon):
        # ...first, old graticule
        loc1 = np.where(lonb1==l)[0] # returns tuple of vectors with dim locations
        if loc1.size==0:
            lon1_idx[i] = lon1_idx[i-1]
        elif loc1.size==1:
            lon1_idx[i] = loc1
        else:
            raise ValueError('Longitude 1 graticule has repeating value %.5f: %s' % (l,lonb1))
        # ...next, new graticule
        loc2 = np.where(lonb2==l)[0] # returns tuple of vectors with dim locations
        if loc2.size==0:
            lon2_idx[i] = lon2_idx[i-1]
        elif loc2.size==1:
            lon2_idx[i] = loc2
        else:
            raise ValueError('Longitude 2 graticule has repeating value %d: %s' % (l,lonb2))
    for i,l in enumerate(lat):
        # ...as above
        loc1 = np.where(latb1==l)[0] # returns tuple of vectors with dim locations
        if loc1.size==0:
            lat1_idx[i] = lat1_idx[i-1]
        elif loc1.size==1:
            lat1_idx[i] = loc1
        else:
            raise ValueError('Latitude 1 graticule has repeating value %d: %s' % (l,latb1))
        # ...next, new graticule
        loc2 = np.where(latb2==l)[0] # returns tuple of vectors with dim locations
        if loc2.size==0:
            lat2_idx[i] = lat2_idx[i-1]
        elif loc2.size==1:
            lat2_idx[i] = loc2
        else:
            raise ValueError('Latitude 2 graticule has repeating value %d: %s' % (l,latb1))
            
    #--------------------------------------------------------------------------
    # Now, use this cell key to downsample (with ugly nested loops)
    #--------------------------------------------------------------------------
    # Get weights -- cosine(lat)*delta_lat for lat, and delta_lon for lon
    # also add newaxis so broadcasting works properly
    latw = (lat[1:]-lat[:-1])*np.abs(np.cos(((lat[1:]+lat[:-1])/2)*np.pi/180))
    lonw = (lon[1:]-lon[:-1])
    latw, lonw = latw[np.newaxis,...]/latw.sum(), lonw[...,np.newaxis]/lonw.sum()
    for i in range(2,data1.ndim):
        # ...need this so that broadcasting rules work properly in loop below
        latw, lonw = latw[...,np.newaxis], lonw[...,np.newaxis]
    # Initialize array
    shape = (lonb2.size-1, latb2.size-1, *data1.shape[2:]) # use tuple expansion here
    data2 = np.empty(shape)
    # Get identifier vectors; for each cell of new graticule, get the indices of
    # old graticule that correspond to each location in combo graticule
    i1 = [lon1_idx[np.where(lon2_idx==i)[0]] for i in range(lonb2.size-1)]
    j1 = [lat1_idx[np.where(lat2_idx==j)[0]] for j in range(latb2.size-1)]
        # so the above gives which square of lons/lats from original graticule correspond 
        # to a given cell in new graticule
    # Get averages
    for i,ii in enumerate(i1):
        for j,jj in enumerate(j1):
            # note because jj, ii are arrays, the indexed array will retain
            # dimensionality even if singleton (recall for a.shape=(3,3), a[1,:]
            # is shape (3,) but a[[1],:] is shape (1,3) [depends if iterable index])
            # also since ii,jj might be different size, must index twice
            data2[i,j,...] = (
                    (data1[ii,:,...][:,jj,...]*latw[:,jj,...]*lonw[ii,:,...]).sum(axis=(0,1)) 
                    /
                    (latw[:,jj,...]*lonw[ii,:,...]).sum()
                    )
                # if numerator is all masked, should match
     
    #--------------------------------------------------------------------------
    # Fix output from previous changes, and return
    #--------------------------------------------------------------------------
    # Remove the duplicate column (was added if new graticule does not have edge on -180/180)
    if flag_cutoff:
        data2 = data2[1:,...]
    # Unroll back to original lon cutoff position
    data2 = np.roll(data2, -lon_roll, axis=0) # unroll
    # Put back into masked array, if necessary
    if flag_masked:
        data2 = np.ma.masked_invalid(data2)
    return data2

def graticule(lon, lat=None, globe=True):
    """
    Gets graticule for full globe of longitude, latitude data.
    Input...
        lon (required): longitude ndarray vector
        lat (required): latitude ndarray vector
        ...or single arg, and do this for arbitrary vector; works in general
        globe (bool, optional):
    ...or
        call with single vector as argument, and function will spit out the halfway
        points plus guesses at the "edges"
    Output...
        latb: ndarray vector of latitude graticule/mesh/grid boundaries
        lonb: ndarray vector of latitude graticule/mesh/grid boundaries
    """
    # enforce monotonic
    if not (lon==np.sort(lon)).all():
        raise ValueError('Longitudes are not monotonic.')
    # fix lons (or generic vector)
    lodiff = lon[1]-lon[0]
    hidiff = lon[-1]-lon[-2]
    lon = np.concatenate(
            (np.array([lon[0]-lodiff/2]), (lon[1:]+lon[:-1])/2, np.array([lon[-1]+hidiff/2]))
            )
    if lat is not None: # this block contains geogrid-specific stuff
        # enforce monotonic
        if not (lat==np.sort(lat)).all():
            raise ValueError('Latitudes are not monotonic.')
        # if grid 'center' reported at pole, publishers really mean from the lower 
        # graticule 'up to' the pole; so modify the reported grid center
        NP = False
        if lat[-1]==90.:
            NP = True
            lat[-1] = 90.-(90.-lat[-2])/4
                # divided by 4, because graticule is halfway between lat[-2] and 90,
                # so we put the new cell center between that line and north pole
        SP = False
        if lat[0]==-90.:
            SP = True
            lat[0] = -90.+(lat[1]+90.)/4
                # same as above
        # fix latitudes
        lodiff = lat[1]-lat[0]
        hidiff = lat[-1]-lat[-2]
        lat = np.concatenate(
                (np.array([lat[0]-lodiff/2]), (lat[1:]+lat[:-1])/2, np.array([lat[-1]+hidiff/2]))
                )
        # and put graticule on NP/SP, if we had that weird ERA-int grid with centers on poles
        if NP or globe: # ...have to enforce because some datasets might go e.g. 84, 86, 88 and STOP; highest cell is 3deg tall
            lat[0] = -90
        if SP or globe:
            lat[-1] = 90
        # warnings
        if (lon[0] % 360) != (lon[-1] % 360) and globe:
            print('Warning: output longitude mesh is not circular.')
    # return stuff
    if lat is None:
        return lon
    else:
        return lon, lat
