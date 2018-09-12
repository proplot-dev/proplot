"""
My netcdf-processing functions; module netCDF4 not necessary, but 
most take as arguments ncgroup variables generated from nc.
"""

#------------------------------------------------------------------------------
# Imports, all
#------------------------------------------------------------------------------
from . import utils
from . import geoutils as geo
import re
import os
import numpy as np
from datetime import datetime
from datetime import timedelta
import netCDF4 as nc4
__all__ = ['ncload', 'ncdump']

#------------------------------------------------------------------------------
# Utilities
#------------------------------------------------------------------------------
def ncload(ncnm, ncgrp_varname, scale=1, add=0,
        time_slice=None,lon_slice=None,lat_slice=None,lev_slice=None, 
        lon_standardize=True, lat_standardize=True, time_standardize=True,
        monthly=False, daily=False, hourly=False,
        lon_pad=False, lon_0=-180, graticule=True,
        squeeze=False, unmask=False, dtype=None,
        verbose=True
        ):
    """
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
    masked: True or False
    
    lon_standardize (bool): standardize lons to be monotonically in range [-180, 180)
    lat_standardize (bool): standardize lats to be monotonically increasing
    time_standardize (bool): get ndarray of datetime objects from time integers, using
        metadata and some fancy string parsing
    monthly (bool): monthly data? set day to 15, truncate hours/minutes/etc.
    daily (bool): daily data? truncate hours/minutes/etc.
    hourly (bool): hourly data? truncate minutes/etc.
    
    unmask (bool): unmask data if we redeived a masked array; False by default,
        but set to True if you know your data has no NaNs
    squeeze (bool): squeeze resulting ndarray/masked array?
    dtype (None or 'type' type): datatype for storing array.

    TODO: Fix .index() problem when input is singleton list.
    TODO: Let ncgrp_varname be a list, and return data as tuple of 
    numpy arrays in that case
    TODO: Fix longitude indexing, so everything is always mod 360... stuff
    can still be messed up it seems.
    """
    
    #--------------------------------------------------------------------------
    # Read file as ncGroup object, and get some info on variable
    #--------------------------------------------------------------------------
    # ...first, the ncgrp
    ncgrp = nc4.Dataset(ncnm, mode='r')
    # ...returns tuple of names (.dimensions) and lengths (.shape)
    filedim_names = ncgrp[ncgrp_varname].dimensions
    filedim_shapes = ncgrp[ncgrp_varname].shape
#     ## Pending, options for multiple varname inputs
#     if type(ncgrp_varname) is list:
#         ncgrp_dimref = ncgrp_varname[0]
#         for iref in ncgrp_varname[1:]:
#             if ncgrp[iref].dimensions != ncgrp[ncgrp_dimref].dimensions:
#                 raise ValueError('Input vars do not have same dimensionality.')
#     else:
#         ncgrp_dimref = ncgrp_varname
#     filedim_names = ncgrp[ncgrp_dimref].dimensions
#     filedim_shapes = ncgrp[ncgrp_dimref].shape

    #--------------------------------------------------------------------------
    # Get dimension locations for "expected" dimension types, also...
    #   1) tools for getting data slice later
    #   2) metadata object with standardized dimension names lon, lat, lev and time
    #--------------------------------------------------------------------------
    # constants
    dim_ndarrayslices = {'lon': lon_slice, 'lat': lat_slice, 'lev': lev_slice, 'time': time_slice}
    dim_ndarraypermute = {'lon':0, 'lat':1, 'lev':2, 'time':3} # fortran/matlab style indexing
    # offer up some 'template' dimension names; see if any exist
    dim_querynames = {
            'lon': ('x','lon','longitude','lo','lons','longs','long','longitudes'),
            'lat': ('y','lat','latitude','la','lats','latitudes'),
            'lev': ('z','lev','level','p','pres','pressure'),
            'time': ('t','time')
            }
    dim_queryresults = {n : [i for i,d in enumerate(filedim_names) if d.lower() in opts]
            for n, opts in dim_querynames.items()} # existence?
    # now get valid locations: where exactly 1 name from each 'template' list is present
    dim_ndarrayids = {n : None if len(val)==0 else -1 if len(val)>1 else val[0]
            for n, val in dim_queryresults.items()} # -1 is a flag
    if any(i==-1 for i in dim_ndarrayids.values()):
        raise ValueError('Ambiguous dimension names.')
    # allow lev to be not found (2D data), but if lon/lev/time not found, raise error
    flag_addlev = False
    dim_badnames = [n for n,v in dim_ndarrayids.items() if v is None] # must construct **list** **beforehand**; can't do inline constructor of any kind, or generator, because we 'pop' an id below
    for n in dim_badnames:
        if n == 'lev':
            dim_ndarrayids.pop('lev')
            dim_ndarrayslices.pop('lev')
            flag_addlev = True
        else:
            raise ValueError('Missing lon, lat, or time dimension.')
    # next get unknown ndarray ids: locations in name tuple that were not filled,
    # and file with the provided FILE name instead of an imposed translation (because we don't know what it is)
    dim_ndarrayids = {**dim_ndarrayids, **{ 
            n : i for i,n in enumerate(filedim_names)
            if i not in dim_ndarrayids.values()
        }}
    # same for slices, use the file names; by default, get all of unknown dimension
    dim_ndarrayslices = {**dim_ndarrayslices, **{
            n : None for i,n in enumerate(filedim_names)
            if i not in dim_ndarrayids.values()
        }} # cat
    # finally, get dictionary for translations
    filedimnames_from_dimnames = { n : filedim_names[dim_ndarrayids[n]] 
            for n in dim_ndarrayids }
    dimnames_from_filedimnames = { n2 : n1 
            for n1, n2 in filedimnames_from_dimnames.items() }
    # and get metadata, with dot notation access
    metadata = utils.DotDict({ dimnames_from_filedimnames[n] : ncgrp.variables[n][:]
            for n in filedim_names })
        # free up for dot notation access!
    
    #--------------------------------------------------------------------------
    #  Try loading time datetimes (use regexps to convert automatically)
    #--------------------------------------------------------------------------
    while time_standardize:
        errmsg = "Failed to convert 'time' to datetime object:"
        # ...get unit string
        if 'unit' in dir(ncgrp[filedimnames_from_dimnames['time']]):
            time_unitstr = ncgrp[filedimnames_from_dimnames['time']].unit
        elif 'units' in dir(ncgrp[filedimnames_from_dimnames['time']]):
            time_unitstr = ncgrp[filedimnames_from_dimnames['time']].units
        else:
            print(errmsg,"ncvar object has no 'unit' or 'units' property.")
            break
        # # pre-process options; remove spaces? bad; e.g. days since 1970-1-1 00:00:00
        # time_unitstr = ' '.join(time_unitstr.split())
        # units (i.e. hours or days; assume days if nothing found)
        if 'minutes' in time_unitstr:
            time_unit = 'minutes'
        elif 'hours' in time_unitstr:
            time_unit = 'hours'
        elif 'days' in time_unitstr:
            time_unit = 'days'
        else:
            print(errmsg,'time unit string (%s) has ambiguous increment unit.' % time_unitstr);
            break
        # find "base" with yyyy-m(m)-d(d) or yyyy/m(m)/d(d) search options; each
        # month and day can have just one digit
        time_fbase = re.search(r'\d{4}[-/]\d{1,2}[-/]\d{1,2}', time_unitstr)
        if time_fbase is not None: 
            time_base = datetime.strptime(
                    time_fbase[0].replace('/','').replace('-',''), '%Y%m%d')
        else:
            print(errmsg,'time unit string (%s) has ambiguous base date.' % time_unitstr);
            break
        time_nums = metadata.time
        # now adjust by converting datenum to timedelta
        if time_unit=='minutes':
            # make vectorized func for elementwise op -- sometimes need to do this 
            # in numpy, for example when another module doesn't play nice with numpy
            time_nums = time_nums.astype(np.int64)*60 # otherwise get SILENT overflow -- scary!
            time_nums = np.vectorize(lambda x: timedelta(seconds=int(x)))(time_nums)
        elif time_unit=='hours':
            # same as above
            time_nums = time_nums.astype(np.int64)*3600 # otherwise get SILENT overflow -- scary!
            time_nums = np.vectorize(lambda x: timedelta(seconds=int(x)))(time_nums)
        elif time_unit=='days':
            # same as above, but arg is days, not seconds
            time_nums = np.vectorize(lambda x: timedelta(days=int(x)))(time_nums)
        if type(time_base) is not datetime:
            time_base = datetime(*time_base)
        # and add timedelta objects to base time
        try:
            time_base + time_nums
        except OverflowError: # sometimes, datasets with null time or monthly climatology
                # will have dateyear '0', which is below datetime.min == (1,1,1,0,0,0)
            time_base = time_base + timedelta(days=366) # since 'year 0' is leap year
            pass
        time = time_base + time_nums # will throw error again if something else is wrong

        # and standardize days for monthly, daily, etc. data
        if monthly:
            time = np.vectorize(lambda x: datetime(x.year, x.month, 15))(time)
        if daily: # eliminate hours/minutes/etc.
            time = np.vectorize(lambda x: datetime(x.year, x.month, x.day))(time)
        if hourly: # eliminate hours/minutes/etc.
            time = np.vectorize(lambda x: datetime(x.year, x.month, x.day, x.hour))(time)
        # and output
        metadata.time = time
        if verbose:
            print(f'Converted time to datetime object:\n' 
                    + f'...original string, {time_unitstr}\n'
                    + f'...base retrieved, {time_base:%B %d, %Y}\n'
                    + f'...unit retrieved, {time_unit}\n'
                    + f'...date range, {metadata.time[0]:%B %d, %Y} to {metadata.time[-1]:%B %d, %Y}')
        break

    #--------------------------------------------------------------------------
    # Deal with longitudes (pre-slicing -- this way, user can reference longitudes
    # from their preferred indexing convention [e.g. 0 to 360] no matter the file convention)
    #--------------------------------------------------------------------------
    if lon_standardize:
        metadata.lon, lon_roll = geo.lonfix(metadata.lon, return_roll=True)
        lon_roll = np.where((metadata.lon % 360 - lon_0 % 360) % 360 == 
                ((metadata.lon % 360 - lon_0 % 360) % 360).min())[0]
        metadata.lon = np.roll(metadata.lon,-lon_roll)
        # ...and fix numbering for monotonic lons
        metadata.lon -= 360*round((metadata.lon[0]-lon_0)/360)
        lon_changeover = np.where(np.diff(metadata.lon)<0)[0]
        if lon_changeover.size!=0:
            metadata.lon[lon_changeover[0]+1:] = (
                    metadata.lon[lon_changeover[0]+1:]+360)
        # ...and temporarily put back
        metadata.lon = np.roll(metadata.lon,lon_roll)

    #--------------------------------------------------------------------------
    # Get dimension ids assciated with slices; apply slices to metadata
    #--------------------------------------------------------------------------
    for n,r in dim_ndarrayslices.items():
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
            if n=='time' and not isinstance(r[0],datetime):
                r[0] = datetime(*r[0])
                if len(r)>1:
                    r[1] = datetime(*r[1])
            # get ids along relevant dimension of dataset ndarray
            if len(r)==1:
                # ...exact location; .where returns tuple with vector
                # for each dimension, so since this is singleton, take first
                dim_ndarrayslices[n] = np.where(getattr(metadata,n)==r[0])[0]
            else:
                # ...between two locations, inclusive
                dim_ndarrayslices[n] = np.where((getattr(metadata,n)>=r[0]) 
                        & (getattr(metadata,n)<=r[1]))[0]
                
            # apply skipping frequency to ids (::n includes endpoints)
            if len(r)==3:
                dim_ndarrayslices[n] = dim_ndarrayslices[n][::r[2]]

            # get actual dimension values associated with ids
            setattr(metadata,n,getattr(metadata,n)[dim_ndarrayslices[n]])

    #--------------------------------------------------------------------------
    # Deal with grid: longitudes (post-slicing) and latitudes
    #--------------------------------------------------------------------------
    # longitudes
    if lon_standardize:
        lon_roll = np.where(metadata.lon==metadata.lon.min())[0]
            # already is ascending in the way we want, from earlier processing
        metadata.lon = np.roll(metadata.lon,-lon_roll)
    if lon_pad:
        metadata.lon = np.pad(metadata.lon,(0,1),'wrap')
        metadata.lon[-1] += 360
    # latitudes; enforce monotonic
    flag_lat = False
    if lat_standardize and len(metadata.lat)>1:
        if metadata.lat[1]<metadata.lat[0]:
            metadata.lat = np.flipud(metadata.lat)
            flag_lat = True
    # built graticule, if requested
    if graticule:
        metadata.lonb, metadata.latb = geo.graticule(metadata.lon, metadata.lat)

    #--------------------------------------------------------------------------
    # Load data with user-requested slices; unmask (if requested)
    #--------------------------------------------------------------------------
    dim_ndarrayslices = {n: range(0,filedim_shapes[dim_ndarrayids[n]]) if dim_ndarrayslices[n]
            is None else dim_ndarrayslices[n] for n in dim_ndarrayids}
    subscripts = [dim_ndarrayslices[dimnames_from_filedimnames[n]] for n in filedim_names]
    data = ncgrp[ncgrp_varname][subscripts]
    # Change type if requested
    if dtype is not None:
        data = data.astype(dtype)
    # Unmask if requested (nc4 loads up masked array if _FillValue or missing_value present)
    if type(data) is np.ma.MaskedArray and unmask:
        data.fill_value = np.nan
        data = data.filled() # does conversion

    #--------------------------------------------------------------------------
    # Finish processing data (add 'z' coord if necessary, premute into imposed
    # dim order, match to lon/lat flipping/slicing/rolling, scale units, squeeze)
    #--------------------------------------------------------------------------
    # fix levels
    if flag_addlev:
        data = data[...,np.newaxis]
        dim_ndarrayids['lev'] = data.ndim-1 # must be -1; that is location in array
        dim_ndarrayslices['lev'] = 1
        metadata.lev = np.array([])
    # permute; append unknown dimensions to the end of array
    data = np.transpose(data, [dim_ndarrayids[n] for n in dim_ndarraypermute]
            + [dim_ndarrayids[n] for n in dim_ndarrayids if n not in dim_ndarraypermute])
    # flipping, circular shifting, and padding
    if flag_lat:
        data = np.flip(data, axis=dim_ndarraypermute['lat'])
    if lon_standardize:
        data = np.roll(data, -lon_roll, axis=dim_ndarraypermute['lon'])
    if lon_pad:
        pad_spec = [(0,0)]*data.ndim
        pad_spec[dim_ndarraypermute['lon']] = (0,1)
        data = np.pad(data, pad_spec, 'wrap')
    # scale data for new units
    data = data*scale + add
    # squeeze, possibly (must do this after processing everything else)
    if squeeze:
        data = data.squeeze()

    #--------------------------------------------------------------------------
    # Return everything
    #--------------------------------------------------------------------------
    if verbose:
        print('Variable "%s" loaded from %s, size %s.' % (ncgrp_varname, ncnm, data.shape))
    return data, metadata

def ncdump(ncnm, verbose=True):
    """
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters...
    ncnm (positional): a netCDF file name, to be loaded by nc4
    verbose (kwarg): Boolean; whether or not nc_attrs, nc_dims, and nc_vars arerinted

    Returns...
    nc_attrs: A Python list of the NetCDF file global attributes
    nc_dims: A Python list of the NetCDF file dimensions
    nc_vars: A Python list of the NetCDF file variables
    """
    # Get file
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
    if verbose:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(ncgrp.getncattr(nc_attr)))
    nc_dims = [dim for dim in ncgrp.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verbose:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim )
            print("\t\tsize:", len(ncgrp.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in ncgrp.variables]  # list of nc variables
    if verbose:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", ncgrp.variables[var].dimensions)
                print("\t\tsize:", ncgrp.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


