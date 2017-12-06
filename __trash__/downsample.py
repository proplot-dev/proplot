def downsample(lonb1, latb1, data1, lonb2, latb2, ignorenan=True):
    import numpy as np
    """
    SHOULD NOT USE ANYMORE; MAKES MORE SENSE TO JUST DOWNSAMPLE THE ORIGINAL NC
    FILE AND SAVE, AS IT IS EXPENSIVE OPERATION. USE NCO FOR THIS.
    ALSO TESTED IT, AND WHILE THIS TAKES E.G. LIKE 10-20 SECONDS TO INTERPOLATE 1BY1 TO 5BY5 (150 YEARS,MONTHLY),
    CDO TAKES 1 SECOND... SO NOT ONLY IS CDO MORE EXTENSIBLE TO OTHER WEIRD GRIDS, BUT IS MUCH FASTER.
    ALSO, FROM .25 RESOLUTION TO 5 RESOLUTION LIKE 20-30 SECONDS, BUT 3 SECONDS FOR CDO
    Downsamples data onto a more sparse grid using area-weighted averagine. User
    should provide the MESHED lonb1/latb1 box boundaries.
    Input...
        lonb1: cell graticule longitudes
        latb1: cell graticule latitudes
        data: data lying on above graticule, shaped lon by lat
        lonb2: interpolant; new cell graticule longitudes
        latb2: interpolant; new cell graticule latitudes
        ignorenan (bool): whether to ignore NaNs; default True (e.g. for SST, want to 
            ignore areas that are 100% land, then in regional average just weight cells by ocean fraction)
    Breakpoint of longitudes DOES NOT MATTER, as long as the graticule is circular
    i.e. the first matches the last mod 360, and lonbitudes are monotonic. Everything
    will be np.roll'ed to match.
    Latitude graticules must also be monotonic.
    """
     
    # Longitudes, latitudes enforce float
    lonb1, latb1 = lonb1.astype(float), latb1.astype(float)
    lonb2, latb2 = lonb2.astype(float), latb2.astype(float)
    if type(data1) is np.ma.MaskedArray:
        flag_masked = True
        data1 = data1.filled(np.nan)
    else:
        flag_masked = False
    # Raise errors
    if data1.shape[0] != lonb1.size-1:
        raise ValueError('Longitude must be on first dimension of data.')
    if data1.shape[1] != latb1.size-1:
        raise ValueError('Latitude must be on second dimension of data.')
    for l in (latb1, latb2):
        if l[0]!=-90. or l[-1]!=90.:
            raise ValueError('Both latitude graticules must start at -90 and '
                    + 'end at 90 latitude.')
    for l in (lonb1, lonb2):
        if (l[0]%360) != (l[-1]%360):
            raise ValueError('Both longitude graticules must be circular.')
     
    # # Fix longitudes to same 360-range, and fill in seams
    # lonb1, data1 = lonfix(lonb1, data1, lonmin=-180)
    # lonmin = lonb2[0] # will permute back later
    # lonb2 = lonfix(lonb2, lonmin=-180)
    #     # will roll back to original position, after downsampling
    # if lonb1[0]>lonmin:
    #     lonb1, data1 = seamfix(lonb1, data1, lonmin=-180)
    # lonb2_old = lonb2.copy()
    # if lonb2[0]>lonmin:
    #     lonb2 = seamfix(lonb2, lonmin=-180)
    #     flag_cutoff = True
    # else:
    #     flag_cutoff = False
    
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
            
    # Get weights -- cosine(lat)*delta_lat for lat, and delta_lon for lon
    # also add newaxis so broadcasting works properly
    latw = (lat[1:]-lat[:-1])*np.abs(np.cos(((lat[1:]+lat[:-1])/2)*np.pi/180))
    lonw = (lon[1:]-lon[:-1])
    latw, lonw = latw[None,...]/latw.sum(), lonw[...,None]/lonw.sum()
    for i in range(2,data1.ndim):
        # ...need this so that broadcasting rules work properly in loop below
        latw, lonw = latw[...,None], lonw[...,None]
    weights = latw*lonw
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
    with np.errstate(divide='ignore', invalid='ignore'): # ignores N/0 and 0/0
        N, percent = len(i1)*len(j1), 0
        print('Downsampling...')
        for i,ii in enumerate(i1): # note fancy indexing with ii/jj retains dimensions
            if ((i+1)*len(j1)/N)>(.01*percent):
                print('%d%% finished...' % (percent,))
                percent = percent+10
            for j,jj in enumerate(j1):
                # ...new
                d = data1[ii,:,...][:,jj,...]
                w = weights[ii,:,...][:,jj,...]
                if ignorenan:
                    f = np.isfinite(d) # filter
                    w = w*f # add to weights
                    d[~f] = 0 # fast assignment
                data2[i,j,...] = (d*w).sum(axis=(0,1))/w.sum(axis=(0,1))
                # ...old
                # data2[i,j,...] = (
                #     (data1[ii,:,...][:,jj,...]*latw[:,jj,...]*lonw[ii,:,...]).sum(axis=(0,1))
                #     /
                #     (latw[:,jj,...]*lonw[ii,:,...]).sum()
                #     )
    data2[np.isinf(data2)] = np.nan # set invalid values to NaN, including /0 stuff
     
    # Final steps; return longitudes to original positions, add back mask
    # if flag_cutoff:
    #     lonb2, data2 = lonb2_old, data2[1:,...]
    # _, data2 = lonfix(lonb2, data2, lonmin=lonmin)
    if flag_masked:
        data2 = np.ma.masked_invalid(data2)
    return data2
