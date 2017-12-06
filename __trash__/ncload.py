def ncload(filename, param,
        times=None, lons=None, lats=None, levs=None, 
        monthoffset=None, # option to manually offset each month; implemented for CESM-LE data
        monthly=False, daily=False,  # at least hourly
        lonmin=-180, # for longitude control
        globe=True, # for graticule
        landfracs=False, # land fractions
        force3d=False, dtype=None, fill_value=None, nanmask=None, # controlling data array
        verbose=True # show what is happening, step by step
        ):
    """
    SHOULD NOT USE ANYMORE; MAKES MORE SENSE TO LOAD INTO XARRAY; IS WELL-DEVELOPED,
    AND UNDERLYING STRUCTURES USES NETCDF4; MUCH MUCH BETTER.

    Loads ndarray of data from .nc file into predictable format and dimension indexing,
    along with metadata structure that can be dot-indexed. Many useful options.

    lons: longitude range, as 2d list/tuple
    lats: latitude range
    levs: level range (for whatever hieght-like index dataset has)
    times: time starting and ending times, as...
        1) length 2 list/tuple of lists/tuples with 3-6 elements
        2) length 2 list/tuple of datetime.datetime objects
      *all slices can have 3rd dimension, indicating a skipping frequency
      for retrieval (useful e.g. for getting 12Z times from 0-6-12-18)
      **ranges will be INCLUSIVE of both endpoints
    
    lonmin: new central longitude; data starts from this minus 180
    globe (bool): globe kwarg passed to my.graticule; will enforce -90 and +90 on 
        ends of lat graticule, even if dataset uses weird indexing.
    
    monthly (bool): monthly DATA? set day to 15, truncate hours/minutes/etc.
    daily (bool): daily DATA? truncate hours/minutes/etc., center to hour 12

    unmask (bool): unmask DATA if we redeived a masked array; False by default,
        but set to True if you know your DATA has no NaNs
    fill_value: *override* fill_value -- add_offset these elements to the mask, if netCDF4
        fails to do so.
    standardize (bool): standardize units to known values, or do not?
    mask: ndarray matching size of data; option to combine with this mask
    force3d (bool): add singleton 'height' dimension data is 2d?
    dtype (None or 'type' type): datatype for storing array.

    TODO: Fix .index() problem when input is singleton list.
    TODO: Let param be a list, and return DATA as tuple of 
    numpy arrays in that case
    TODO: Fix longitude indexing, so everything is always mod 360... stuff
    can still be messed up it seems.
    """
#     # Pending, options for multiple varname inputs
#     if type(param) is list:
#         ncgrp_dimref = param[0]
#         for iref in param[1:]:
#             if file[iref].dimensions != file[ncgrp_dimref].dimensions:
#                 raise ValueError('Input vars do not have same dimensionality.')
#     else:
#         ncgrp_dimref = param
#     dim_filenames = file[ncgrp_dimref].dimensions
#     dim_fileshapes = file[ncgrp_dimref].shape
     
    # THINGS COMPLETELY AUTOMATED BY XARRAY LOAD_DATASET
    # Read the file
    M = Dict() # my special dictionary
    with nc4.Dataset(filename, mode='r') as file:
        # Get variable info
        # ...returns tuple of names (.dimensions) and lengths (.shape)
        NCparam  = file[param]
        ncnames  = NCparam.dimensions
        ncshapes = NCparam.shape
        # ...retrieve units, name
        M.name  = NCparam.long_name
        M.units = NCparam.units

        # Get dimension locations for "expected" dimension types, also...
        slices = dict(lon=lons, lat=lats, lev=levs, time=times)
        positions = ('lon', 'lat', 'lev', 'time')
        # ...offer up some 'template' dimension names; see if any exist
        querynames = dict(time=('t','time', 'date'),
                lon = ('x','i','lon','longitude','lo','lons','longs','long','longitudes'), # have seen i/j in files
                lat = ('y','j','lat','latitude','la','lats','latitudes'),
                lev = ('z','lev','level','p','pres','pressure'), # t is 'transport-z'
                )
        queryresults = {n: [i for i,d in enumerate(ncnames) if d.lower() in opts]
                for n, opts in querynames.items()} # existence?
        # ...now get valid locations: where exactly 1 name from each 'template' list is present
        # ...then raise errors for duplicate names, names not found
        ncpositions = {n: None if len(val)==0 else -1 if len(val)>1 else val[0]
                for n, val in queryresults.items()} # -1 is a flag
        if any(i==-1 for i in ncpositions.values()): raise ValueError('Ambiguous dimension names.')
        is3d, hastime = True, True # by default, assume 3d data with time
        badnames = [n for n,v in ncpositions.items() if v is None] # unknown
        for n in badnames:
            if n == 'lev':
                is3d = False
                slices.pop('lev')
                ncpositions.pop('lev')
            elif n=='time':
                hastime = False
                slices.pop('time')
                ncpositions.pop('time')
            else:
                raise ValueError('Missing x or y dimension.')
        # ...next, add to dictionaries the UNKNOWN dimension names (using the provided name); slicing
        # the unknown dimension names is not allowed
        ncpositions.update({n:i for i,n in enumerate(ncnames) if i not in ncpositions.values()})
        slices.update({n:slice(None) for i,n in enumerate(ncnames) if i not in ncpositions.values()})
        
        # Translate from the NetCDF dimension names to my standardized names, and vice versa
        mynames2ncnames = {n: ncnames[ncpositions[n]] for n in ncpositions}
        ncnames2mynames = {n2: n1 for n1, n2 in mynames2ncnames.items()}
        
        # Add variables to metadata dictionary
        M.update({ncnames2mynames[n]: file[n][:] for n in ncnames})
        
        # Interpret time dimension info
        if hastime:
            # Get time unit string, and calendar
            NCtime = file[mynames2ncnames['time']]
            units = NCtime.units
            try: calendar = NCtime.calendar # e.g. noleap, for no leap years
            except AttributeError: calendar = 'gregorian' # default; currently used

            # Use netCDF4.num2date to convert to datetimes, extremely awesome function that even allows for special calendar types
            M.time = nc4.num2date(M.time, units, calendar=calendar) # calendar may be gregorian, julian, 'noleap', etc.
                # this creates ndarray of datetime objects if Gregorian calendar, otherwise
                # creates ndarray of pseudo-datetime object that have same behavior, with extra attributes dayofwk, dayofyr
            try: M.time = np.array([t._to_real_datetime() for t in M.time])
            except AttributeError: pass # when already gives us datetimes
                # converts back to normal datetime object; we lose the dayofwk, dayofyr attributes

            # Standardize datetimes for daily/monthly/etc.; sometimes monthly data might have day on last day of month, etc.
            if monthly:
                M.time = np.array([datetime(t.year, t.month, 15, 0) for t in M.time]) # centered on month
            elif daily:
                M.time = np.array([datetime(t.year, t.month, t.day, 12) for t in M.time]) # centered on day
            else: # by default, data is atmospheric data is at least hourly
                M.time = np.array([datetime(t.year, t.month, t.day, t.hour) for t in M.time]) # point obs

            # Special treatment for stupid dumb terribly organized datasets (CESM-LE)
            if monthoffset is not None:
                pseudomonths = np.array([t.month + monthoffset for t in M.time]) # for -1 offset, get 0-11 for example
                M.time = np.array([datetime(t.year + (p-1)//12, ((p-1) % 12) + 1, t.day, t.hour) for t,p in zip(M.time, pseudomonths)])
                    # error will likely be raised if this is not monthly data; e.g. February 31st
                    # note that this removes the dayofwk, dayofyr information
            if verbose: print('Datetime string "%s", retrieved %s to %s.' % 
                    (units, M.time[0].strftime('%B %d %Y'), M.time[-1].strftime('%B %d %Y')))

        # Slicing: Get dimension ids assciated with slices; apply slices to M
        M.lon = geo.lonfix(M.lon, lonmin=lonmin, roll=False) # make sure lons span expected 360-degree range
        ncslices = {n:slice(None) for n in slices} # initialize slice object dict
        for n,s in slices.items():
            if s is not None:
                try: iter(s) # a location
                except TypeError: s = (s, s) # a range/range+step
                if len(s) not in (2,3):
                    raise ValueError('Slice must be length 1 (point), 2 (range), or 3 (range, with skipping interval).')
                loc = np.where((M[n]>=s[0]) & (M[n]<=s[1]))[0] # inclusive region
                if len(s)==2:
                    ncslices[n] = slice(loc.min(),loc.max()+1) # no skipping
                elif len(s)==3:
                    ncslices[n] = slice(loc.min(), loc.max(), s[2]) # skipping
                M[n] = M[n][ncslices[n]] # slice up dimensions vectors

        # Deal with grid: latitudes and longitudes
        # ...longitudes (get new grid longitude minimum after slicing)
        lon_roll = -np.argmin(M.lon) # the 360-range convention is already set
        M.lon = np.roll(M.lon, lon_roll)
        # ...latitudes; enforce monotonic
        flag_lat = False
        if M.lat[1]<M.lat[0]:
            M.lat, flag_lat = np.flipud(M.lat), True

        # Retrieve DATA with user-requested slices
        # ...load, using slices; just unpack dictionary and order according to dimension order in filename
        DATA = NCparam[[ncslices[ncnames2mynames[n]] for n in ncnames]]

    # Apply permuting/rolling stuff, 
    # unmask (if requested), and process masked array properties
    # ...convert from mask (keep it simple, don't need that)
    if type(DATA) is np.ma.MaskedArray:
        DATA = DATA.filled(np.nan)
        if verbose: print(f'Dataset unmasked.')
    # ...hange type if requested (must be done before anything else)
    if dtype is not None:
        DATA = DATA.astype(dtype)
        if verbose: print(f'Datatype set to {dtype}.')
    # ...manually override fill value
    if fill_value is not None:
        DATA[DATA.round()==fill_value] = np.nan
        if verbose: print(f'Manually masked data matching integer fill_value {fill_value}.')
    # ...fix levels
    if not is3d:
        DATA = DATA[...,None] # new level dimension is on end of axes
        ncpositions['lev'] = DATA.ndim-1 # must be -1; that is location in array
        M.lev = np.array([])
    if not hastime:
        DATA = DATA[...,None] # new level dimension is on end of axes
        ncpositions['time'] = DATA.ndim-1 # must be -1; that is location in array
        M.time = np.array([])
    # ...permute; append unknown dimensions to the end of array
    DATA = np.transpose(DATA, [ncpositions[n] for n in positions] # order of positions == order of sort
            + [ncpositions[n] for n in ncpositions if n not in positions])
    if not is3d and not force3d:
        DATA = DATA.squeeze()
    # ...flipping, circular shifting, and padding
    if flag_lat: DATA = np.flip(DATA, axis=positions.index('lat'))
    DATA = np.roll(DATA, lon_roll, axis=positions.index('lon'))
    if verbose: print('Longitudes rolled; new leftmost cell is closest above %d.' % lonmin)
    # ...add mask, if requested
    if nanmask is not None:
        # ...get np.tile tuple (e.g. mask over ocean, and data has time dimension at end)
        mask_tile = (1,)*nanmask.ndim + DATA.shape[nanmask.ndim:]
        for i in range(DATA.ndim-nanmask.ndim):
            nanmask = nanmask[...,None]
        nanmask = np.tile(nanmask, mask_tile)
        # ...then apply new nanmask
        DATA[nanmask] = np.nan
        if verbose: print('Custom NaN mask applied.')

    # Manage data units, and other extra stuff
    # ...load land fractions, if requested or required by another option
    if landfracs: # ...must do this after closing file
        print('Loading land fractions at grid resolution...')
        M.landfracs = geo.landfracs(M.lon, M.lat)
    if verbose: print('Data units: "%s".' % M.units)

    # Return everything
    if verbose: print('Parameter "%s" loaded from %s, size %s.' % (param, filename, DATA.shape))
    return DATA, M

