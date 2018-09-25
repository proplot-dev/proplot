## My netcdf-processing functions; module netCDF4 not necessary, but 
## most take as arguments ncgroup variables generated from nc$

def ncdate(datenums, unit='days', base=(1900,1,1,0,0)):
    '''
    Returns a datetime object np array; datenums should be np.array format.
    Handles everything without dealing with division, decimals, etc. using datetime 
    and timedelta classes, so very very safe!
    '''
    import numpy as np
    from datetime import datetime
    from datetime import timedelta
    # convert to ndarray of timedelta objects
    # note a good way to convert datetime->timedelta is dt-datetime.min;
    # equivalent of matlab datenum, but in more neat timedelta object
    if unit=='hours':
        # make vectorized func for elementwise op -- sometimes need to do this 
        # in numpy, for example when another module doesn't play nice with numpy
        datenums = datenums.astype(np.int64)*3600 # otherwise get SILENT overflow -- scary!
        datenums = np.vectorize(lambda x: timedelta(seconds=int(x)))(datenums)
    elif unit=='days':
        # same as above, but arg is days, not seconds
        datenums = np.vectorize(lambda x: timedelta(days=int(x)))(datenums)
    else:
        raise ValueError("Unknown unit: " + str(unit) + "Use 'hours' or 'days'.")
    # now get actual dates
    if type(base) is not datetime:
        base = datetime(*base)
    return base + datenums
    # ...simple as that! numpy just puts base datetime into array, then trys
    # an elementwise plus operator on everything, and it works!

def ncload(ncnm, ncgrp_varname, scale=1, add=0,
        time_slice=None,lon_slice=None,lat_slice=None,lev_slice=None, 
        lon_standardize=True, lat_standardize=True, time_standardize=True,
        lon_pad=False, lon_0=-180):
    '''
    Sets up the data for plotting; accounts for units, etc.  Will enforce
    monotonic longitude from the specied lon_0.

    scale: multiplicative scaling from default units, to desired units
    add: additive scaling from default units (comes after *scale)

    lon_slice: longitude range, as 2d list/tuple
    lat_slice: latitude range
    lev_slice: level range (for whatever hieght-like index dataset has)
    time_slice: time starting and ending times, as...
        1) length 2 list/tuple of lists/tuples with 3-6 elements
        2) length 2 list/tuple of datetime.datetime objects
      *all slices can have 3rd dimension, indicating a skipping frequency
      for retrieval (useful e.g. for getting 12Z times from 0-6-12-18)
      **ranges will be INCLUSIVE of both endpoints

    lon_pad: pad longitudes around breaking point, for smooth contours on
        certain map projections (e.g. polar azimuthal)
    lon_0: cutoff longitude; 1st in array is lowest longitude >= this value

    TODO: Fix .index() problem when input is singleton list.
    TODO: Let ncgrp_varname be a list, and return data as tuple of 
    numpy arrays in that case
    TODO: Fix longitude indexing, so everything is always mod 360... stuff
    can still be messed up it seems.
    '''
    ## Imports
    import re
    import numpy as np
    from datetime import datetime
    import netCDF4 as nc4
    
    ## Read file
    ncgrp = nc4.Dataset(ncnm, mode='r')

#     ## Pending, options for multiple varname inputs
#     if type(ncgrp_varname) is list:
#         ncgrp_dimref = ncgrp_varname[0]
#         for iref in ncgrp_varname[1:]:
#             if ncgrp[iref].dimensions != ncgrp[ncgrp_dimref].dimensions:
#                 raise ValueError('Input vars do not have same dimensionality.')
#     else:
#         ncgrp_dimref = ncgrp_varname
#     ncgrp_vardims = ncgrp[ncgrp_dimref].dimensions
#     ncgrp_varshape = ncgrp[ncgrp_dimref].shape
    ## Load dimension stuff
    ncgrp_vardims = ncgrp[ncgrp_varname].dimensions
    ncgrp_varshape = ncgrp[ncgrp_varname].shape

    ## Initial stuff
    #dim_order = ['lon','lat','lev','time'] # matlab style
    dim_order = ['time','lev','lat','lon'] # python uses c-indexing; get used to it
    dim_slices = {'lon': lon_slice, 'lat': lat_slice, 'lev': lev_slice, 'time': time_slice}
    dim_testnms = {
            'lon': ['x','lon','longitude','lo','lons','longs','long','longitudes'],
            'lat': ['y','lat','latitude','la','lats','latitudes'],
            'lev': ['z','lev','level','p','pres','pressure'],
            'time': ['t','time']
            }

    ## Get dimension locations for "expected" dimension types
    dim_ids = {n : [i for i,d in enumerate(ncgrp_vardims) if d.lower() in opts]
            for n, opts in dim_testnms.items()} 
    dim_ids = {n : None if len(val)==0 else -1 if len(val)>1 else val[0]
            for n, val in dim_ids.items()} # -1 is a flag
    if any(i==-1 for i in dim_ids.values()):
        raise ValueError('Ambiguous dimension names.')

    ## Unknown and missing dimensions; consider, and modify dim_ids, dim_slices
    dim_missing = [n for n,v in dim_ids.items() if v is None]
    flag_addlev = False
    for n in dim_missing:
        if n == 'lev':
            dim_ids.pop('lev')
            dim_slices.pop('lev')
            flag_addlev = True
        else:
            raise ValueError('Missing lon, lat, or time dimension.')
    dim_idsunknown = {n : i for i,n in enumerate(ncgrp_vardims)
            if i not in dim_ids.values()}
    dim_slicesunknown = {n : None for i,n in enumerate(ncgrp_vardims)
            if i not in dim_ids.values()}
    dim_ids = {**dim_ids, **dim_idsunknown} # cat
    dim_slices = {**dim_slices, **dim_slicesunknown} # cat
        # note we use our "imposed" names for expected dimensions, but use
        # the default ones for unknown dimensions
    dim_truenms = { n : ncgrp_vardims[dim_ids[n]] for n in dim_ids }
    dim_truenms_inverse = { n2 : n1 for n1, n2 in dim_truenms.items() }

    ## Load metadata
    metadata = { n : ncgrp.variables[n][:] for n in ncgrp_vardims }
    for newnm,oldnm in dim_truenms.items():
        metadata[newnm] = metadata.pop(oldnm)
            # pop returns old contents upon request

    ## Deal with times (use regexps to convert automatically)
    while time_standardize:
        errmsg = "Failed to convert 'time' to datetime object."
        # ...get unit string
        if 'unit' in dir(ncgrp[dim_truenms['time']]):
            time_unitstr = ncgrp[dim_truenms['time']].unit
        elif 'units' in dir(ncgrp[dim_truenms['time']]):
            time_unitstr = ncgrp[dim_truenms['time']].units
        else:
            print(errmsg);
            break
        # units (i.e. hours or days; assume days if nothing found)
        if 'hours' in time_unitstr:
            time_unit = 'hours'
        elif 'days' in 'time_unitstr':
            time_unit = 'days'
        else:
            print(errmsg);
            break
        # find "base", with two search options: yyyy/mm/dd and yyyy-mm-dd
        # looks gross, but cleanest way actually and extensible for new date formats
        time_qbase = [r'\d{4}-\d{2}-\d{2}', r'\d{4}/\d{2}/\d{2}']
        time_fbase = [s for s in (re.search(s,time_unitstr) for s in time_qbase)
                if s is not None] # failure returns None; the parentheses makes a temporary iterable
        # retrieve base from singleton list containing regexp object
        if time_fbase: # list non-empty
            time_base = datetime.strptime(
                    time_fbase[0][0].replace('/','').replace('-',''), '%Y%m%d')
            metadata['time'] = ncdate(
                    metadata['time'], unit=time_unit, base=time_base)
            break
        else:
            print(errmsg);
            break

    ## Deal with longitudes
    if lon_standardize:
        lon_roll = np.where((metadata['lon'] % 360 - lon_0 % 360) % 360 == 
                ((metadata['lon'] % 360 - lon_0 % 360) % 360).min())[0]
        metadata['lon'] = np.roll(metadata['lon'],-lon_roll)
        # ...and fix numbering for monotonic lons
        metadata['lon'] -= 360*round((metadata['lon'][0]-lon_0)/360)
        lon_changeover = np.where(np.diff(metadata['lon'])<0)[0]
        if lon_changeover.size!=0:
            metadata['lon'][lon_changeover[0]+1:] = (
                    metadata['lon'][lon_changeover[0]+1:]+360)
        # ...and temporarily put back
        metadata['lon'] = np.roll(metadata['lon'],lon_roll)

    ## Process slices
    for n,r in dim_slices.items():
        if r is not None:
            # if not iterable (e.g. single number), place into singleton list
            try:
                iter(r)
            except:
                r = [r] # e.g., just a number
                pass
            r = list(r) # need item assignment sometimes
            # test length of dimension slice input
            if not 1 <= len(r) <= 3:
                raise ValueError('Slice must be length 1 (point) or '
                    + '2 (range), or 3 (range, with skipping interval).')
            # enforce equivalent types of starting/ending locations
            if len(r)>1:
                if type(r[0]) is not type(r[1]):
                    raise ValueError('Slice range types are inconsistent.')
            # convert time list/tuple ends to datetime
            if n == 'time' and not isinstance(r[0],datetime):
                r[0] = datetime(*r[0])
                if len(r)>1:
                    r[1] = datetime(*r[1])
            # get ids along relevant dimension of dataset ndarray
            if len(r)==1:
#                 # old method... what are you smoking?
#                 dim_slices[n] = metadata[n].index(r)
                # ...exact location; .where returns tuple with vector
                # for each dimension, so since this is singleton, take first
                dim_slices[n] = np.where(metadata[n]==r[0])[0]
            else:
#                 # old method... what are you smoking?
#                 dim_slices[n] = [list(metadata[n]).index(v) for v in
#                         list(metadata[n]) if r[0] <= v <= r[1]]
                # ...between two locations, inclusive
                dim_slices[n] = np.where((metadata[n]>=r[0]) 
                        & (metadata[n]<=r[1]))[0]
                
            # apply skipping frequency to ids (::n includes endpoints)
            if len(r)==3:
                dim_slices[n] = dim_slices[n][::r[2]]
#             # old method... what are you smoking? unnecessary list creation
#             metadata[n] = [m for i,m in enumerate(metadata[n])
#                     if i in dim_slices[n]]
            # get actual dimension values associated with ids
            metadata[n] = metadata[n][dim_slices[n]]

    ## Deal with longitudes
#     metadata['lon'] = np.array(metadata['lon']) # re-cast as array, if filtered
    if lon_standardize:
        lon_roll = np.where(metadata['lon']==metadata['lon'].min())[0]
            # already is ascending in the way we want, from earlier processing
        metadata['lon'] = np.roll(metadata['lon'],-lon_roll)
    if lon_pad:
        metadata['lon'] = np.pad(metadata['lon'],(0,1),'wrap')
        metadata['lon'][-1] += 360

    ## Deal with latitudes (enfore monotonic)
    flag_lat = False
#     metadata['lat'] = np.array(metadata['lat']) # re-cast as array, if filtered
    if lat_standardize and len(metadata['lat'])>1:
        if metadata['lat'][1]<metadata['lat'][0]:
            metadata['lat'] = np.flipud(metadata['lat'])
            flag_lat = True

    ## Load data with user-requested slices
    dim_slices = {n: range(0,ncgrp_varshape[dim_ids[n]]) if dim_slices[n]
            is None else dim_slices[n] for n in dim_ids}
    subscripts = [dim_slices[dim_truenms_inverse[n]] for n in ncgrp_vardims]
    data = ncgrp[ncgrp_varname][subscripts]
#    if all(d is None for d in dim_slices):
#        data = ncgrp[ncgrp_varname][:] # just take everything

    ## Scale data for new units
    data = data*scale + add

    ## Fix levels
    if flag_addlev:
        data = data[...,np.newaxis]
        dim_ids['lev'] = data.ndim
        dim_slices['lev'] = 1
        metadata['lev'] = {}

    ## Permute into imposed dimension order
    data = np.transpose(data, [dim_ids[n] for n in dim_order]
            + [dim_ids[n] for n in dim_ids if n not in dim_order])

    ## Fix data for lon, lat
    if flag_lat:
        # for some reason, flip doesn't exist
        data = data.swapaxes(0,dim_order.index('lat'))
        data = data[::-1,...]
        data = data.swapaxes(0,dim_order.index('lat'))
    if lon_standardize:
        data = data.swapaxes(0,dim_order.index('lon'))
        data = np.roll(data,-lon_roll,axis=0)
        data = data.swapaxes(0,dim_order.index('lon'))
        #data = np.roll(data,-lon_roll,axis=4)
    if lon_pad:
        data = data.swapaxes(0,dim_order.index('lon'))
        data = np.pad(data,[(1,0)] + [(0,0)]*(data.ndim-1),'wrap')
        data = data.swapaxes(0,dim_order.index('lon'))
        #pad_spec = [(0,0)]*data.ndim
        #pad_spec[dim_order.index('lon')] = (1,0)
        #data = np.pad(data,pad_spec,'wrap')

    ## Wrap things up
    print('Variable %s loaded.' % ncgrp_varname)
    return data, metadata

def ncdump(ncnm, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    ncgrp : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars arerinted

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    # Imports
    import netCDF4 as nc4
    import os
    print(os.getcwd())
    print(os.listdir())
    ncgrp = nc4.Dataset(ncnm, mode='r')
    
    # Definition
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:",repr(ncgrp.variables[key].dtype))
            for ncattr in ncgrp.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,repr(ncgrp.variables[key].getncattr(ncattr)))
        except KeyError:
           print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = ncgrp.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(ncgrp.getncattr(nc_attr)))
    nc_dims = [dim for dim in ncgrp.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim )
            print("\t\tsize:", len(ncgrp.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in ncgrp.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", ncgrp.variables[var].dimensions)
                print("\t\tsize:", ncgrp.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


