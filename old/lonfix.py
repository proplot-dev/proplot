def lonfix(lon, data=None, lonmin=-180, roll=True):
    """
    FALLS INTO EXAMPLE OF SOMETHING THAT DOESN'T DESERVE ITS OWN FUNCTION; REALLY 
    ONLY NEED TO CONSIDER THIS WHEN DOWNLOADING, AND WHEN PLOTTING.
    Shifts longitudes to the breakpoint specified by lonmin.
     * roll (bool): whether we roll to lonmin, or just enforce monotonicity from lonmin.
     * lonmin: the new minimum longitude; lon[0] <-- lon[np.argmin(lon>=lonmin-180)]
     * return_roll: whether to return the rolling parameter
    """
    # Initial stuff, and checks
    lon = np.array(lon) # enforce array
    lon_old = lon.copy()
    if lon.max()>lon.min()+360:
        raise ValueError('Longitudes must span 360 degrees at most.')
    if lon.min()<-360 or lon.max()>360:
        raise ValueError('Longitudes must fall in range [-360, 360].')
    if lonmin<-360 or lonmin>0:
        raise ValueError('Minimum longitude must fall in range [-360, 0].')
    
    # Enforce new longitude convention
    lon -= 720
    while True:
        filter_ = lon<lonmin
        if filter_.sum()==0:
            break
        lon[filter_] += 360
    if roll:
        lon_roll = -np.argmin(lon) # always returns FIRST value; if go from 0...359,0 (borders), will
            # return initial value
        if lon[0]==lon[-1]:
            lon = np.roll(lon[:-1], lon_roll)
            lon = np.append(lon, lon[0]+360)
        else:
            lon = np.roll(lon, lon_roll)
        if data is not None:
            data = np.roll(data, lon_roll, axis=0)

    # Return stuff, and remember rules of tuple expansion
    if data is None:
        return lon
    else:
        return lon, data

