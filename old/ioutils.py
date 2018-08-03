"""
My netcdf-processing functions; module netCDF4 not necessary, but 
most take as arguments ncgroup variables generated from nc.
"""

#------------------------------------------------------------------------------
# Imports and global functions
#------------------------------------------------------------------------------
# from .basics import Dict
from glob import glob
import xarray as xr
import dateutil.parser as dparser # date parser
import regex as re # includes additional functionality
import os
import subprocess, threading
import ecmwfapi as ecmwf
import numpy as np
from datetime import date, datetime, timedelta
import calendar
import netCDF4 as nc4
# __all__ = ['eraint', 'ncload', 'ncdump']
DATADIR = os.environ['HOME'] + '/data/'
LANDFRACS = DATADIR + '/landfracs.nc'
    # for land fractions, and other shapefiles that might be added
    # to this function

#------------------------------------------------------------------------------
# Helper classes
#------------------------------------------------------------------------------
class Dict(dict):
    """
    Dot notation access to dictionary attributes; EXTREMELY useful, and no
    slower than normal dictionaries.
    """
    # Very simple
    # __getattr__ = dict.get
    # __setattr__ = dict.__setitem__
    # __delattr__ = dict.__delitem__
    # More complex implimentation
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v
        if kwargs: # if non-empty
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr) # make getattr() call dict().get()

    def __setattr__(self, key, value):
        self.__setitem__(key, value) # make setattr() call dict().__setitem__()

    def __setitem__(self, key, value):
        super().__setitem__(key, value) # make simple assignment set the attr, and update the dict
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item) # make delattr use dictionary delete

    def __delitem__(self, key):
        super().__delitem__(key) # make del command delete dictionary item
        del self.__dict__[key]

class Metadata():
    def __init__(self):
        # Class for assigning metadata
        # Must create this... or forget it, because xarray is way better
        pass

#------------------------------------------------------------------------------
# Utilities
#------------------------------------------------------------------------------
def eraint(params, stream, levtype,
        daterange=None, yearrange=None, monthrange=None, dayrange=None,
        years=None, months=None, # can specify list
        levrange=None, levs=None,
        hours=(0,6,12,18),
        res=1.0, box=None,
        filename='eraint.nc'):
    """
    Retrieves ERA-Interim DATA using the provided API. User MUST have, in home
    directory, a file named '.ecmwfapirc'; see API documentation, but should look like...
        {
        "url"   : "https://api.ecmwf.int/v1",
        "key"   : "960dbe61271d3902c8b0f768d69d679f",
        "email" : "email@gmail.com"
        }
    ...with the key found on your user/profile page on the ecmwf website.
    
    Time range arguments:
    years/yearrange -- list of range of years
    months/monthrange -- list or range of months
    daterange -- range of dates/datetimes

    Other input arguments:
    params: can be either of...
        list/tuple of variable string names
        individual variable string
    *** Must know MARS id for requested params; can add to dictionary in code below using
        https://rda.ucar.edu/datasets/ds627.0/docs/era_interim_grib_table.html
    stream: can be any of...
        'synoptic' (6 hourly DATA)
        'monthly' (monthly mean)
    levtype: can be any of...
        'pl' (pressure levels)
        'sfc' (earth surface)
        'pt' (potential temperature)
        'pv' (2pvu surface)
    levrange: can be either of...
        length-2 tuple/list of pressure/pt levels; retrieves all available levels between these
        single number, to pick individual level
    levs: can be either of...
        list/tuple of multiple levels; retrieves each level in list
        single number, to pick individual level 
    hours: can be either of...
        list/tuple of integer hours (should be in [0,6,12,18])
        single number, to pick individual hour
    res: desired output resolution; not sure if this is arbitrary, or if ERA-interim only has
        a select few valid resolution options.
    box: can be either of...
        string name for particular region, e.g. "europe" (see documentation)
        the West/South/East/North boundaries (so lower-left corner, upper-right corner), as a length-4 list/tuple
    filename: name of file output
    """
    # Data stream
    stream = { # oper is original, moda is monthly mean of daily means
            'synoptic':  'oper',
            'monthly':   'moda'
            }.get(stream)
    if stream is None:
        raise ValueError('Must choose from "oper" for synoptic fields, "moda" for monthly means of daily means.')

    # Variable id conversion (see: https://rda.ucar.edu/datasets/ds627.0/docs/era_interim_grib_table.html)
    if isinstance(params, str):
        params = (params,)
    params = [{
            't2m':     '167.128', # 2m temp
            'd2m':     '168.128', # 2m dew point
            'sst':     '34.128', # sst
            'msl':     '151.128', # sea level pressure
            'slp':     '151.128', # same
            'z':       '129.128', # geopotential
            't':       '130.128', # temp
            'u':       '131.128', # u wind
            'v':       '132.128', # v wind
            'w':       '135.128', # w wind
            'q':       '133.128', # specific humidity
            'r':       '157.128', # relative humidity
            'vort':    '138.128', # relative vorticity
            'vo':      '138.128', # same
            'zeta':    '138.128', # same
            'pt':      '3.128', # potential temp (available on 2pvu surf)
            'theta':   '3.128', # same
            'p':       '54.128', # pressure (availble on pt, 2pvu surfaces)
            'pres':    '54.128', # same
            'pv':      '60.128', # potential vorticity (available on p, pt surfaces)
            'precip':  '228.128',
            }.get(p) for p in params] # returns generator object for each param
    if None in params:
        raise ValueError('MARS id for variable is unknown (might need to be added to this script).')
    params = '/'.join(params)
    
    # Time selection as various RANGES or LISTS
    # ...priority; just use daterange as datetime or date objects
    if daterange is not None:
        try: iter(daterange)
        except TypeError:
            daterange = (daterange,) # want a SINGLE DAY
        # options for monthly or daily data
        if stream=='moda':
            y0, m0, y1, m1 = daterange[0].year, daterange[0].month, daterange[1].year, daterange[1].month
            N = max(y1-y0-1, 0)*12 + (13-m0) + m1 # number of months in range
            dates = '/'.join('%04d%02d00' % (y0 + (m0+n-1)//12, (m0+n-1)%12 + 1) for m in range(N))
        else:
            dates = '/to/'.join(d.strftime('%Y%m%d') for d in daterange) # MARS will get calendar days in range
    # ...alternative; list the years/months desired, and if synoptic, get all calendar days within
    else:
        # ...first, years
        if years is not None:
            try: iter(years)
            except TypeError:
                years = (years,) # single month
        elif yearrange is not None:
            try: iter(yearrange)
            except TypeError: # single year
                years = (yearrange,)
            else:
                years = tuple(range(yearrange[0], yearrange[1]+1))
        else:
            raise ValueError('You must use "years" or "yearrange" kwargs.')
        # ...next, months (this way, can just download JJA data, for example)
        if months is not None:
            try: iter(months)
            except TypeError:
                months = (months,) # single month
        elif monthrange is not None:
            try: iter(monthrange)
            except TypeError: # single year
                months = (monthrange, monthrange)
            else:
                months = tuple(range(monthrange[0], monthrange[1]+1))
        else:
            months = tuple(range(1,13))
        # ...and get dates; options for monthly means and daily stuff
        if stream=='moda':
            dates = '/'.join(
                '/'.join('%04d%02d00' % (y,m) for m in months)
                for y in years)
        else:
            dates = '/'.join(
                '/'.join(
                '/'.join('%04d%02d%02d' % (y,m,i+1) for i in range(calendar.monthrange(y,m)[1]))
                for m in months)
                for y in years)
            
    # Level selection as RANGE or LIST
    # ...update this list if you modify script for ERA5, etc.
    levchoices = {
            'sfc':  None,
            'pv':   None,
            'pl':   np.array([1,2,3,5,7,10,20,30,50,70,100,125,150,175,200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000]),
            'pt':   np.array([265,270,285,300,315,330,350,370,395,430,475,530,600,700,850]),
            }.get(levtype, [])
    if levchoices==[]:
        raise ValueError('Invalid level type. Choose from "pl", "pt", "pv", "sfc".')
    if levtype not in ('sfc','pv'): # these have multiple options
        # require input
        if levs is None and levrange is None and levtype not in ('sfc','pv'):
            raise ValueError('Must specify list of levels to "levs" kwarg, range of levels to "levrange" kwarg, or single level to either one.')
        # convert levels to mars request
        if levs is not None:
            try: iter(levs)
            except TypeError: # single level
                levs = (levs,)
        elif levrange is not None:
            try: iter(levrange)
            except TypeError: # single level
                levs = (levrange,)
            else:
                levs = levchoices[(levchoices>=levrange[0]) & (levchoices<=levrange[1])].flat
        levs = '/'.join(str(l) for l in levs)

    # Other parameters
    # Resolution
    res = '%.5f/%.5f' % (res,res) # same in latitude/longitude required, for now
    # Area - can be specified as pre-defined region (e.g. string 'europe') OR n/s/w/e boundary
    if box is not None and type(box) is not str:
        box = '/'.join(str(b) for b in (box[3], box[0], box[2], box[1]))
    # Hour conversion
    try: iter(hours)
    except TypeError:
        hours = (hours,)
    hours = '/'.join(str(h).zfill(2) for h in hours) # zfill padds 0s on left

    # Server instructions
    # ...not really sure what happens in some situations: list so far...
    # 1) evidently if you provide with variable string-name instead of numeric ID, 
    #       MARS will search for correct one; if there is name ambiguity/conflict will throw error
    # 2) on GUI framework, ECMWF only offers a few resolution options, but program seems
    #       to run when requesting custom resolutions like 5deg/5deg
    retrieve = {
        'class':    'ei', # ecmwf classifiction; choose ERA-Interim
        'expver':   '1',
        'dataset':  'interim', # ...thought we already did that; *shrug*
        'type':     'an', # type of field; analysis 'an' or forecast 'fc'
        'resol':    'av', # prevents truncation before transformation to geo grid
        'step':     '0', # number of hours forecast has been run into future from 'time'
        'gaussian': 'reduced',
        # 'format':   'netcdf',
        # 'grid':     res,
        'format':   'grib', # can also spit raw output into GRIB; apparently
            # ERA-Interim uses bilinear interpolation to make grid of point obs,
            # which makes sense, because their reanalysis model just picks out point observations
            # from spherical harmonics... so maybe grid cell concept is dumb? maybe need to focus
            # on just using cosine weightings, forget about rest?
        # 'grid':     'N32',
        'stream':   stream, # product monthly, raw, etc.
        'date':     dates,
        'time':     hours,
        'levtype':  levtype,
        'param':    params,
        'target':   filename, # save location
        }
    if levs is not None: retrieve.update(levelist=levs)
    if box is not None: retrieve.update(area=box)
    if stream!='moda': retrieve.update(hour=hour)
    print('Final MARS request: %s' % retrieve)

    # Retrieve DATA with settings
    server = ecmwf.ECMWFDataServer()
    server.retrieve(retrieve)
    return

# def cdo(shfile, *args):
def cdo(cdo, *args, **kwargs):
    """
    Calls CDO process, with arguments specified by args.
    Just a convenience function so I can call CDO methods from python scripts.
    Will use the builtin python-cdo module, which has careful exception handling, and
    in my testing it seemed just using subprocess sometimes resulted in errors for commands
    that worked from shell; parses out the 'process' and args.
    This looks much nicer in practical use than using the Cdo module directly; its 
    implementation for long chains is really ugly.
    Examples:
        sub (subtracting climate)
        ensmean (ensemble mean, many files)
        ensaverage (ensemble mean, ignoring nan)
        enspctl,50 (ensemble median, many files; use --percentile=numpy to interpolate)
        -select,name=<name> (for selecting individual variable from file)
        -sellonlatbox,W,E,S,N (for selecting a subset region)
            note only operators beginning with DASH can be chained; the rest
            have variable number of input files, so cannot; although, can
            put the first one as file
    NOTE infiles can also be list of PIPES/FILES; operator chaining with -.
    """
    # Extract options
    options = ''
    while len(args)>0:
        if not args[0].startswith('-'):
            break
        else:
            options += args[0] # add string
            args = args[1:]
        if len(args)==0:
            raise ValueError('Invalid input.')
    # Parse input
    if type(cdo) is str:
        raise ValueError('First argument must be an initialized Cdo() instance.')
    print('CDO command: ' + ' '.join(args))
    arg1, input_ = args[0].split(','), args[1:]
    process, args = arg1[0], arg1[1:] # if args was singleton, args[1:] is empty

    # Call cdo command
    output = getattr(cdo, process)(*args, 
            input=input_, options=options,
            **kwargs # output options; includes None, returnArray, returnCdf
            )

    # # Just call CDO command, and raise error
    # print('CDO command: ' + ' '.join(args))
    # subprocess.check_call(( 'cdo', '-C', '-O', # colorized, overwrite
    #     # '-v', # verbose
    #     '-L', '--no_warnings', # suppress warnings, and lock input/output; otherwise can get segfaults for big chains
    #     *args), close_fds=True) # arguments
    
    # And return outfile (the last argument, always)
    if 'output' in kwargs:
        print('Result sent to output file: %s' % kwargs['output'])
    return output

def mergedname(infiles, fmt='%Y%m'):
    """
    Assigns name to file composed of other merged files.
    """
    # List files
    if isinstance(infiles, str):
        infiles = glob(infiles)
        # want to look at the filenames explicitly, rename based on them
    
    # Isolate date part of each file
    regex = (fmt.replace('%Y','\d{4}').replace('%m','\d{2}').replace('%d','\d{2}') # date
            .replace('%H','\d{2}').replace('%M','\d{2}').replace('%M','\d{2}')) # time
    if '%' in regex:
        raise ValueError('Unrecognized format string. Use only numeric Y/m/d/H/M/S identifiers.')
    counts = [len(re.findall(regex, f, overlapped=True)) for f in infiles]
                # remember generators destroy stuff only after you look at particular item
    if any(c>1 for c in counts):
        print('\n'.join(infiles))
        raise IOError('Found too many matches, or time-merged netCDF file already exists, in at least one of listed files.')
    if any(c==0 for c in counts):
        print('\n'.join(infiles))
        raise IOError('Could not find match for at least one of listed files.')
    
    # Figure out outfile name
    times = [datetime.strptime(re.search(regex, f).group(), fmt) for f in infiles]
    datestring = min(times).strftime(fmt) + '-' + max(times).strftime(fmt)
    outfile = infiles[0].replace(re.search(regex, infiles[0]).group(), datestring)
        # use first filename as template
    
    # And exit, with message
    return outfile

def gridfile(gridfile, lon, lat):
    """
    Creates formatted file suitable for input as "gridfile" into a CDO interpolation
    scheme, from the vectors of longitude/latitude centers.
    """
    # Check input
    # loninc, latinc = lon[1]-lon[0], lat[1]-lat[0]
    loninc, latinc = np.diff(lon).mean(), np.diff(lat).mean()
        # might have tiny floating point differences...
    
    # Run CDO method
    with open(gridfile, 'w') as grid:
            # the w+ makes new file, if does not exist
            # also w+ overwrites existing content; a+ only appends
        # may need to overwrite lon/lat names
        grid.write('\n'.join((
            'gridtype = lonlat',
            'gridsize = %d' % (lon.size*lat.size),
            'xsize    = %d' % lon.size,
            'ysize    = %d' % lat.size,
            'xfirst   = %.5f' % lon[0], # gives .1 arcsecond precision
            'yfirst   = %.5f' % lat[0],
            'xinc     = %.5f' % loninc,
            'yinc     = %.5f' % latinc,
            # 'xname  = lon',
            # 'yname  = lat',
            ))) # combines strings with '\n' newlines
        
    # Return
    return gridfile

def nc2xr(filename, param, lonmin=-180,
        times=None, lons=None, lats=None, levs=None, # for slicing, specify coordinate/tuple range/length-3 tuple range + skip
        engine='netcdf4', cache=False, # don't load into ndarray, until explicitly asked
        ):
    """
    Loads xarray object from netCDF data.
    """
    # Load dataset
    with xr.open_dataset(filename) as file:
        # Time should already be fixed.
        # ...do nothing
        
        # Standardize remaining dimension names (to avoid making copies, must
        # be done on Dataset object, not Dataarray object; the latter has no inplace option)
        querynames = dict( # time must ALWAYS be named 'time'
                lon = ('x','i','lon','longitude','lo','long'), # have seen i/j in files
                lat = ('y','j','lat','latitude','la'),
                lev = ('z','lev','level','p','pres','pressure','pt','theta','h','height'), # t is 'transport-z'
                ) # standardized enough, that this should be pretty safe
        for assignee, options in querynames.items():
            count = 0
            for opt in options:
                if opt in file.indexes:
                    file.rename({opt: assignee}, inplace=True)
                    count += 1
            if count==0:
                if assignee not in ('lev','time'):
                    raise IOError('Candidate for "%s" dimension not found.' % assignee)
                # else:
                #     print('Warning: data has no "%s" dimension.' % assignee)
                #         # allow this, sometimes
            elif count>1:
                raise IOError('Multiple candidates for "%s" dimension found.' % assignee)
            
        # Get data
        data = file[param]

        # Slice up, with sel
        def makeslice(arg):
            if arg is None: 
                arg = (None,)
            try: iter(arg)
            except TypeError:
                if arg is None: arg = slice(None)
            else:
                arg = slice(*arg)
            return arg # a scalar, Noneslice (select all, or ':'), or range
                # don't use singleton slices of scalars, because these silently
                # fail and select all whereas a single slice would raise KeyError if index not present
        slices = dict(
                lon = makeslice(lons),
                lat = makeslice(lats)
                )
        if 'lev' in data.indexes:
            slices.update(lev = makeslice(levs))
        if 'time' in data.indexes:
            slices.update(time = makeslice(times))
        data = data.sel(**slices) # slice dat shit up yo
            # this action also copies it from filesystem, it seems
        
        # Finally, load view of DataArray on disk into memory
        data.load()
    
    # Fix precision of time units; sometimes get weird useless round-off error; 
    # then restore to numpy datetime64[ns] because xarray seems to require it
    # ALSO ran into mysterious problem where dataset could be loaded, but then
    # COULD NOT BE SAVED because one of the datetimes wasn't serializable... this was
    # in normal data, the CCSM4 CMIP5 results, made no sense; range was 1850-2006
    if 'time' in data.indexes:
        data['time'] = data.time.values.astype('datetime64[D]').astype('datetime64[ns]')

    # Make latitudes monotonic (note extracting values way way faster)
    if data.lat.values[0]>data.lat.values[1]:
        data = data.isel(lat=slice(None,None,-1))

    # Enforce longitude ordering convention
    # Not always necessary, but this is safe/fast; might as well standardize
    values = data.lon.values - 720 # equal mod 360
    while True:
        filter_ = values<lonmin
        if filter_.sum()==0:
            roll = values.argmin()
            data = data.roll(lon=-roll)
            data['lon'] = np.roll(values, -roll)
            break
        values[filter_] += 360

    # Re-order dims to my expected order
    # order = ['lon', 'lat', 'lev', 'time']
    order = ['time','lev','lat','lon'] # better for numpy broadcasting
    if 'lev' not in data.indexes:
        order.remove('lev')
    if 'time' not in data.indexes:
        order.remove('time')
    data.transpose(*order)
    
    # Wanted to do all of the above with 
    return data # will remove singleton dimensions

def nc2cube(filename, param,
        times=None, lons=None, lats=None, levs=None, # for slicing, specify coordinate/tuple range/length-3 tuple range + skip
        landfracs=False, # don't load into ndarray, until explicitly asked
        ):
    """
    Loads netCDF data into iris cube; think about implementing this.
    """
    # WRITE CODE
    # MORE CODE
    # HERE HAVE SOME CODE
    return

def ncload(filename, param,
        times=None, lons=None, lats=None, levs=None, 
        monthoffset=None, # option to manually offset each month; implemented for CESM-LE data
        monthly=False, daily=False,  # at least hourly
        lonmin=-180, # for longitude control
        landfracs=False, # land fractions
        force3d=False, dtype=None, fill_value=None, nanmask=None, # controlling data array
        verbose=True # show what is happening, step by step
        ):
    """
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
        M.lon -= 720
        while True:
            filter_ = M.lon<lonmin
            if filter_.sum()==0: break
            M.lon[filter_] += 360
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
        if lon_roll!=0: M.lon = np.roll(M.lon, lon_roll)
        # ...latitudes; enforce monotonic
        flag_lat = False
        if M.lat[1]<M.lat[0]: M.lat, flag_lat = np.flipud(M.lat), True

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
    # if not is3d and not force3d:
    #     DATA = DATA.squeeze()
    # ...flipping, circular shifting, and padding
    if flag_lat: DATA = np.flip(DATA, axis=positions.index('lat'))
    if lon_roll!=0: DATA = np.roll(DATA, lon_roll, axis=positions.index('lon'))
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
        f = remapcon(LANDFRACS, M.lon, M.lat) # file on disk
        M.landfracs = ncload(f, 'data')[0].squeeze()
    if verbose: print('Data units: "%s".' % M.units)

    # Return everything
    if verbose: print('Parameter "%s" loaded from %s, size %s.' % (param, filename, DATA.shape))
    return DATA, M

def ncdump(filename, verbose=True):
    """
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters...
    filename (positional): a netCDF file name, to be loaded by nc4
    verbose (kwarg): Boolean; whether or not nc_attrs, nc_dims, and nc_vars arerinted

    Returns...
    nc_attrs: A Python list of the NetCDF file global attributes
    nc_dims: A Python list of the NetCDF file dimensions
    nc_vars: A Python list of the NetCDF file variables
    """
    # Definition
    def print_attr(key):
        try:
            print("\t\ttype:", file.variables[key].dtype)
            for ncattr in file.variables[key].ncattrs():
                print('\t\t%s:' % ncattr, file.variables[key].getncattr(ncattr))
        except KeyError:
           print("\t\tWARNING: %s does not contain variable attributes" % key)

    # Load global attributes, dim info, var info (inside with statement, so it closes afterward)
    with nc4.Dataset(filename, mode='r') as file:
        nc_attrs = file.ncattrs()
        nc_dims = [dim for dim in file.dimensions]  # list of nc dimensions
        nc_vars = [var for var in file.variables]  # list of nc variables
        # Print info
        if verbose:
            print("NetCDF Global Attributes:")
            for nc_attr in nc_attrs:
                print('\t%s:' % nc_attr, file.getncattr(nc_attr))
            print("NetCDF dimension information:")
            for dim in nc_dims:
                print("\tName:", dim )
                print("\t\tsize:", len(file.dimensions[dim]))
                print_attr(dim)
            print("NetCDF variable information:")
            for var in nc_vars:
                if var not in nc_dims:
                    print('\tName:', var)
                    print("\t\tdimensions:", file.variables[var].dimensions)
                    print("\t\tsize:", file.variables[var].size)
                    print_attr(var)
    return nc_attrs, nc_dims, nc_vars
