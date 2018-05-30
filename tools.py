#------------------------------------------------------------------------------
# All
#------------------------------------------------------------------------------
import os
import numpy as np
try:
    import ecmwfapi as ecmwf
except ModuleNotFoundError:
    print("Warning: ECMWF API unavailable.")
import xarray as xr
import scipy.signal as signal
import scipy.stats as stats
import time as T # have some big algorithms here; want to guage how long they take
from . import const
# from . import geoutils as geo
# from . import basics
# __all__ = [
#     'deriv1', 'deriv2', 'deriv_uneven', 'diff', 'laplacian',
#     'absvo', 'pt', 'pv', 'qgpv', 'eqlat', 'waq', 'waqlocal',
#     'power', 'cps', 'lowpass', 'lanczos',
#     ]

#------------------------------------------------------------------------------
# Helper functions, for permuting/flattening/unflattening/unpermuting things
#------------------------------------------------------------------------------
def _flatten(data, end=2):
    """
    Flatten trailing dimensions onto single dimension (useful).
    Call this to manipulate data on first dimensions, then restore with unflatten.
    * Note that numpy prod of an empty iterable will be 1; so, adds singleton dim.
    * Note we use tuple expansion of the [shape] tuple
    """
    return np.reshape(data, (*data.shape[:end], np.prod(data.shape[end:]).astype(int)), order='F'), data.shape

def _unflatten(data, shape, end=2):
    """
    Undo action of flatten; determine new shape with tuple expansion.
    """
    if data.shape[-1] != np.prod(shape[end:]):
        raise ValueError('Number of trailing dimensions %d does not match trailing shape %s.' % (data.shape[-1], shape[end:]))
    return np.reshape(data, (*data.shape[:end], *shape[end:]), order='F')

def _permute(data, axis=-1):
    """
    Permutes a given axis onto the LAST dimension.
    """
    return np.rollaxis(data, axis, data.ndim)

def _unpermute(data, axis=-1, select=-1):
    """
    Undoes action of permute; permutes LAST dimension back onto original axis.
    Rolls axis <select> until it lies "before" (i.e. to the left) of axis <axis>.
    """
    if select<0: select = data.ndim+select
    return np.rollaxis(data, select, axis)

#-------------------------------------------------------------------------------
# NetCDF loading tools
#-------------------------------------------------------------------------------
def nc(filename, param, lonmin=0, times=None, lons=None, lats=None, levs=None, **kwargs):
    """
    Function for loading up NetCDF data, and standardizing dimension names and order,
    and standardizing longitude range (i.e. -180-180 or 0-360). Extremely simple.
    * Also converts time dimension from np.datetime64 to list of python datetime.datetime objects.
    * For slicing, specify coordinate/tuple range/length-3 tuple range + skip
    """
    # Helper function; create slice
    def makeslice(arg):
        if arg is None:
            arg = (None,)
        try: iter(arg)
        except TypeError:
            pass
        else:
            arg = slice(*arg)
        return arg # a scalar, Noneslice (select all, or ':'), or range
    # Load dataset
    if not os.path.exists(filename):
        raise ValueError(f'{filename} does not exist!')
    with xr.open_dataset(filename, engine='netcdf4', cache=False, **kwargs) as file:
        # Standardize remaining dimension names (to avoid making copies, must
        # be done on Dataset object, not Dataarray object; the latter has no inplace option)
        querynames = { # time must ALWAYS be named 'time'
                'lat':  ['y','j','lat','latitude','la'],
                'lon':  ['x','i','lon','longitude','lo','long'], # have seen i/j in files
                'lev':  ['z','lev','level','p','pres','pressure','pt','theta','h','height'], # t is 'transport-z'
                } # standardized enough, that this should be pretty safe
        for assignee, options in querynames.items():
            count = 0
            for opt in options:
                if opt in file.indexes:
                    file.rename({opt: assignee}, inplace=True)
                    count += 1
            if count==0:
                if assignee not in ('lev','time'):
                    raise IOError(f'Candidate for "{assignee}" dimension not found.')
            elif count>1:
                raise IOError(f'Multiple candidates for "{assignee}" dimension found.')
        slices = { 'lon': makeslice(lons),
                   'lat': makeslice(lats) }
        data = file[param]
        if 'lev' in data.indexes:
            slices.update(lev=makeslice(levs))
        if 'time' in data.indexes:
            slices.update(time=makeslice(times))
        data = data.sel(**slices) # slice dat shit up yo
            # this action also copies it from filesystem, it seems
        data.load() # load view of DataArray from disk into memory

    # Fix precision of time units... some notes:
    # 1) sometimes get weird useless round-off error; convert to days, then restore to numpy datetime64[ns]
    #   because xarray seems to require it
    # 2) ran into mysterious problem where dataset could be loaded, but then
    #   COULD NOT BE SAVED because one of the datetimes wasn't serializable... this was
    #   in normal data, the CCSM4 CMIP5 results, made no sense; range was 1850-2006
    if 'time' in data.indexes:
        data['time'] = data.time.values.astype('datetime64[D]').astype('datetime64[ns]')
    # Enforce longitude ordering convention
    # Not always necessary, but this is safe/fast; might as well standardize
    values = data.lon.values-720 # equal mod 360
    while True: # loop only adds 360s to longitudes
        filter_ = values<lonmin
        if filter_.sum()==0: # once finished, write new longitudes and roll
            roll = values.argmin()
            data = data.roll(lon=-roll)
            data['lon'] = np.roll(values, -roll)
            break
        values[filter_] += 360
    # Make latitudes monotonic (note extracting values way way faster)
    try: data.lat.values[1]
    except IndexError:
        pass
    else:
        if data.lat.values[0]>data.lat.values[1]:
            data = data.isel(lat=slice(None,None,-1))

    # Re-order dims to my expected order before returning
    order = ['time','lev','lat','lon'] # better for numpy broadcasting
    if 'lev' not in data.indexes:
        order.remove('lev')
    if 'time' not in data.indexes:
        order.remove('time')
    data.transpose(*order)
    return data

#-------------------------------------------------------------------------------
# ERA-interim
#-------------------------------------------------------------------------------
def eraint(params, stream, levtype,
        daterange=None, yearrange=None, monthrange=None, dayrange=None,
        years=None, months=None, # can specify list
        levrange=None, levs=None,
        hours=(0,6,12,18),
        res=1.0, box=None,
        filename='eraint.nc'):
    """
    Retrieves ERA-Interim DATA using the provided API. User MUST have, in home
    directory, a file named '.ecmwfapirc'; see API documentation, but should look like:
        {
        "url"   : "https://api.ecmwf.int/v1",
        "key"   : "960dbe61271d3902c8b0f768d69d679f",
        "email" : "email@gmail.com"
        }
    with the key found on your user/profile page on the ecmwf website.

    Time range arguments:
    years/yearrange -- list of range of years
    months/monthrange -- list or range of months
    daterange -- range of dates/datetimes

    Other input arguments:
    params: can be either of:
        list/tuple of variable string names
        individual variable string
    *** Must know MARS id for requested params; can add to dictionary in code below using
        https://rda.ucar.edu/datasets/ds627.0/docs/era_interim_grib_table.html
    stream: can be any of:
        'synoptic' (6 hourly DATA)
        'monthly' (monthly mean)
    levtype: can be any of:
        'pl' (pressure levels)
        'sfc' (earth surface)
        'pt' (potential temperature)
        'pv' (2pvu surface)
    levrange: can be either of:
        length-2 tuple/list of pressure/pt levels; retrieves all available levels between these
        single number, to pick individual level
    levs: can be either of:
        list/tuple of multiple levels; retrieves each level in list
        single number, to pick individual level
    hours: can be either of:
        list/tuple of integer hours (should be in [0,6,12,18])
        single number, to pick individual hour
    res: desired output resolution; not sure if this is arbitrary, or if ERA-interim only has
        a select few valid resolution options.
    box: can be either of:
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
    # Priority; just use daterange as datetime or date objects
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
    # Alternative; list the years/months desired, and if synoptic, get all calendar days within
    else:
        # First, years
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
        # Next, months (this way, can just download JJA data, for example)
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
        # And get dates; options for monthly means and daily stuff
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
    # Update this list if you modify script for ERA5, etc.
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
    # Not really sure what happens in some situations: list so far:
    # 1) evidently if you provide with variable string-name instead of numeric ID,
    #       MARS will search for correct one; if there is name ambiguity/conflict will throw error
    # 2) on GUI framework, ECMWF only offers a few resolution options, but program seems
    #       to run when requesting custom resolutions like 5deg/5deg
    retrieve = {
        'class':    'ei', # ecmwf classifiction; choose ERA-Interim
        'expver':   '1',
        'dataset':  'interim', # thought we already did that; *shrug*
        'type':     'an', # type of field; analysis 'an' or forecast 'fc'
        'resol':    'av', # prevents truncation before transformation to geo grid
        'step':     '0', # number of hours forecast has been run into future from 'time'
        'gaussian': 'reduced',
        # 'format':   'netcdf',
        # 'grid':     res,
        'format':   'grib', # can also spit raw output into GRIB; apparently
            # ERA-Interim uses bilinear interpolation to make grid of point obs,
            # which makes sense, because their reanalysis model just picks out point observations
            # from spherical harmonics; so maybe grid cell concept is dumb? maybe need to focus
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

#------------------------------------------------------------------------------
# Miscellaneous/Simple
#------------------------------------------------------------------------------
def year(dt):
    """
    Gets year from numpy datetime object (used e.g. by xarray, pandas).
    """
    return dt.astype('datetime64[Y]').astype(np.int32)+1970 # the astype(int) is actually super fast (ns)
        # and the above in general is faster than list comprehension with container of datetime objects
        # UNIX time starts at 1970-01-01 00:00:00

def month(dt):
    """
    Gets month from numpy datetime object (used e.g. by xarray, pandas).
    """
    return dt.astype('datetime64[M]').astype(np.int32)%12 + 1
        # will convert datetime64 units from [ns] (default) to months, then spit out months relative to year
        # UNIX time starts at 1970-01-01 00:00:00

def arange(min_, *args):
    """
    Duplicate behavior of np.arange, except with inclusive endpoints; dtype is
    controlled very carefully, so should be 'most precise' among min/max/step args.
    Input:
        stop
        start, stop, [step]
        just like np.arange
    Output:
        the array sequence
    """
    # Optional arguments just like np.arange
    if len(args)==0:
        max_ = min_
        min_ = 0 # this re-assignes the NAME "min_" to 0
        step = 1
    elif len(args)==1:
        max_ = args[0]
        step = 1
    elif len(args)==2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer? Get new "max"
    if min_//1==min_ and max_//1==max_ and step//1==step:
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += 1
    # Input is float or mixed; cast all to float64, then get new "max"
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += step/2
        # max_ = np.nextafter(max_, np.finfo(np.dtype(np.float64)).max)
            # gives the next FLOATING POINT, in direction of the second argument
            # forget this; round-off errors from continually adding step to min mess this up
    return np.arange(min_, max_, step)

def match(*args):
    """
    Match arbitrary number of 1D vectors; will return slices for producing the matching
    segment from either vector, and the vector itself, so use as follows:
        i1, i2, ..., vmatch = match(v1, v2, ...)
        v1[i1] == v2[i2] == ... == vmatch
    Useful e.g. for matching the time dimensions of 3D or 4D variables collected
    over different years and months.
    """
    vs = [np.array(v) for v in args]
    if not all(np.all(v==np.sort(v)) for v in vs):
        raise ValueError('Vectors must be sorted.')
    # Get common minima/maxima
    min_all, max_all = max(v.min() for v in vs), min(v.max() for v in vs)
    try:
        min_locs = [np.where(v==min_all)[0][0] for v in vs]
        max_locs = [np.where(v==max_all)[0][0] for v in vs]
    except IndexError:
        raise ValueError('Vectors do not have matching maxima/minima.')
    slices = [slice(min_i, max_i+1) for min_i,max_i in zip(min_locs,max_locs)]
    if any(v[slice_i].size != vs[0][slices[0]].size for v,slice_i in zip(vs,slices)):
        raise ValueError('Vectors are not identical between matching minima/maxima.')
    elif any(not np.all(v[slice_i]==vs[0][slices[0]]) for v,slice_i in zip(vs,slices)):
        raise ValueError('Vectors are not identical between matching minima/maxima.')
    return slices + [vs[0][slices[0]]]
# def match(v1, v2):
#     """
#     Match two 1D vectors; will return slices for producing the matching
#     segment from either vector, and the vector itself, so use as follows:
#         i1, i2, vmatch = match(v1, v2)
#         v1[i1] == v2[i2] == vmatch
#     Useful e.g. for matching the time dimensions of 3D or 4D variables collected
#     over different years and months.
#     """
#     v1, v2 = np.array(v1), np.array(v2)
#     if not np.all(v1==np.sort(v1)) or not np.all(v2==np.sort(v2)):
#         raise ValueError('Vectors must be sorted.')
#     # Get common minima/maxima
#     min12, max12 = max(v1.min(), v2.min()), min(v1.max(), v2.max())
#     try:
#         min1f, min2f = np.where(v1==min12)[0][0], np.where(v2==min12)[0][0]
#         max1f, max2f = np.where(v1==max12)[0][0], np.where(v2==max12)[0][0]
#     except IndexError:
#         raise ValueError('Vectors do not have matching maxima/minima.')
#     slice1, slice2 = slice(min1f, max1f+1), slice(min2f, max2f+1)
#     if v1[slice1].size != v2[slice2].size:
#         raise ValueError('Vectors are not identical between matching minima/maxima.')
#     elif not (v1[slice1]==v2[slice2]).all():
#         raise ValueError('Vectors are not identical between matching minima/maxima.')
#     return slice1, slice2, v1[slice1]

#------------------------------------------------------------------------------
# Geographic
#------------------------------------------------------------------------------
# Handy functions in degrees
def sin(x):
    return np.sin(x*np.pi/180)
def cos(x):
    return np.cos(x*np.pi/180)
def tan(x):
    return np.tan(x*np.pi/180)
def arcsin(x):
    return np.arcsin(x)*180/np.pi
def arccos(x):
    return np.arccos(x)*180/np.pi
def arctan(x):
    return np.arctan(x)*180/np.pi

# Other stuff
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
        # Switch
        has90 = True if lat[-1]==90 else False
        hasm90 = True if lat[0]==-90 else False # need these switches later
        # Use corrections for dumb grids with 'centers' at poles
        if hasm90:
            self.latc[0], self.latb[0] = -90+dlat1/4, -90
        if has90:
            self.latc[-1], self.latb[-1] = 90-dlat2/4, 90
        # Corrected grid widths (cells half as tall near pole)
        self.dlon = self.lonb[1:]-self.lonb[:-1]
        self.dlat = self.latb[1:]-self.latb[:-1]
        if hasm90:
            self.dlat[0] /= 2
        if has90:
            self.dlat[-1] /= 2

        # Theta/phi coordinates
        self.phic, self.phib, self.dphi = self.latc*np.pi/180, self.latb*np.pi/180, self.dlat*np.pi/180
        self.thetac, self.thetab, self.dtheta = self.lonc*np.pi/180, self.lonb*np.pi/180, self.dlon*np.pi/180

        # Area weights (function of latitude only)
        # Includes the latitude correction
        self.weights = self.dphi[None,:]*self.dtheta[:,None]*np.cos(self.phic[None,:])
        self.areas = self.dphi[None,:]*self.dtheta[:,None]*np.cos(self.phic[None,:])*(const.a**2)
        # areas = dphi*dtheta*np.cos(phic)*(const.a**2)[None,:]
            # close approximation to area; cosine is extremely accurate
            # make lon by lat, so broadcasting rules can apply

    def __repr__(self): # what happens when viewing in interactive session
        # if __str__ not set (the print() result), will use __repr__
        n = 3
        return "Properties of grid latitudes/longitudes.\n"\
        f"Lat centers (latc): {', '.join(f'{self.latc[i]:.2f}' for i in range(n))}, ... {self.latc[-1]:.2f}\n"\
        f"Lat borders (latb): {', '.join(f'{self.latb[i]:.2f}' for i in range(n))}, ... {self.latb[-1]:.2f}\n"\
        f"Lat widths (dlat): {', '.join(f'{self.dlat[i]:.2f}' for i in range(n))}, ... {self.dlat[-1]:.2f}\n"\
        f"Lon centers (lonc): {', '.join(f'{self.lonc[i]:.2f}' for i in range(n))}, ... {self.lonc[-1]:.2f}\n"\
        f"Lon borders (lonb): {', '.join(f'{self.lonb[i]:.2f}' for i in range(n))}, ... {self.lonb[-1]:.2f}\n"\
        f"Lon widths (dlon): {', '.join(f'{self.dlon[i]:.2f}' for i in range(n))}, ... {self.dlon[-1]:.2f}\n"\
        "Coordinates in radians also provided (lat=phi, lon=theta).\n"\
        "Approximate grid cell areas also provided (longitude x latitude)."

def geopad(lon, lat, data, nlon=1, nlat=0):
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

def geomean(lon, lat, data, box=(None,None,None,None),
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
    a = Properties(lon, lat).areas

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

    Input:
    lon1, lat1, lon2, lat2 (positional): ndarrays of longitude and latitude in degrees
    * Each **pair** should have identical shape
    * If both pairs are non-scalar, they must **also** be identically shaped
    * If one pair is scalar, distances are calculated between the scalar pair
        and each element of the non-scalar pari.
    Output:
    d: ndarray of distances in km
    """
    # Earth radius, in km
    R = const.a*1e-3
    # Convert to radians, get differences
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    # Haversine, in km
    km = 2.*R*np.arcsin(np.sqrt(
        np.sin(dlat/2.)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.)**2
        ))
    return km

#------------------------------------------------------------------------------
# Statistics
#------------------------------------------------------------------------------
def gaussian(mean, sigma):
    """
    Returns sample points on Gaussian curve.
    """
    norm = stats.norm(loc=mean, scale=sigma)
    x = np.linspace(norm.ppf(0.0001), norm.ppf(0.9999), 1000) # get x through percentile range
    pdf = norm.pdf(x)
    return x, pdf

#------------------------------------------------------------------------------
# Math/simple
#------------------------------------------------------------------------------
def rolling(data, window=5, axis=-1):
    """
    Read this: https://stackoverflow.com/a/4947453/4970632
    Generates rolling numpy window along final axis; can then operate with
    functions like polyfit or mean along the new last axis of output.
    Just creates *view* of original array, without duplicating data, so no worries
    about efficiency.
    * Will generate a new axis in the -1 position that is a running representation
      of value in axis numver <axis>.
    * Strides are apparently the 'number of bytes' one has to skip in memory
      to move to next position *on the given axis*. For example, a 5 by 5
      array of 64bit (8byte) values will have array.strides == (40,8).
    * Should consider using swapaxes instead of these _permute and _unpermute
      functions, might be simpler.
    TODO: Add option that preserves *edges* i.e. does not reduces length
    of dimension to be 'rolled' by (window-1).
    """
    # Roll axis, reshape, and get generate running dimension
    # data = np.rollaxis(data, axis, data.ndim)
    if axis<0: axis = data.ndims+axis # e.g. if 3 dims, and want to axis dim -1, this is dim number 2
    data = _permute(data, axis)
    shape = data.shape[:-1] + (data.shape[-1]-(window-1), window)
    strides = [*data.strides, data.strides[-1]] # repeat striding on end
    data = np.lib.stride_tricks.as_strided(data, shape=shape, strides=strides)
    return _unpermute(data, axis, select=-2) # want to 'put back' axis -2;
        # axis -1 is our new "rolling"/"averaging" dimension, keep that in same position

def running(*args, **kwargs):
    """
    Easier to remember name.
    """
    return rolling(*args, **kwargs)

def roots(poly):
    """
    Find real-valued root for polynomial with input coefficients; pretty simple.
    Format of input array p: p[0]*x^n + p[1]*x^n-1 + ... + p[n-1]*x + p[n]
    """
    # Just use numpy's roots function, and filter results.
    r = np.roots(poly) # input polynomial; returns ndarray
    filt = (r==np.real(r))
    r = r[filt].astype(np.float32) # take only real-valued ones
    return r

def slope(x, y, axis=-1):
    """
    Get linear regression along axis, ignoring NaNs. Uses np.polyfit.
    """
    # First, reshape
    y = np.rollaxis(y, axis, y.ndim).T # reverses dimension order
    yshape = y.shape # save shape
    y = np.reshape(y, (yshape[0], np.prod(yshape[1:])), order='F')
    # Next, supply to polyfit
    coeff = np.polyfit(x, y, deg=1) # is ok with 2d input data
        # DOES NOT accept more than 2d; also, much faster than loop with stats.linregress
    # And finally, make
    slopes = coeff[0,:]
        # Coefficients are returned in reverse; e.g. 2-deg fit gives c[0]*x^2 + c[1]*x + c[2]
    slopes = np.reshape(slopes, (1, y.shape[1:]), order='F').T
    return np.rollaxis(slopes, y.ndim-1, axis)

#------------------------------------------------------------------------------
# Finite-difference schemes
#------------------------------------------------------------------------------
# def derivl(h, y, axis=0, accuracy=2, keepedges=True):
#     """
#     Differentiation using left-assignment.
#     Will keep edges by default (this is supposed to be quick and easy method).
#     Should forget this and just use centered differencing instead.
#     """
#     # Simple Euler scheme
#     # y = np.rollaxis(y, axis, y.ndim)
#     y = _permute(y, axis)
#     diff = (y[...,1:]-y[...,:-1])/h
#     if keepedges:
#         diff = np.concat((diff, diff[...,-1:]), axis=-1) # repeat it
#     # return np.rollaxis(diff, y.ndim-1, axis)
#     return _unpermute(diff, axis)

def deriv1(h, y, axis=0, accuracy=2, keepleft=False, keepright=False, keepedges=False):
    """
    First order finite differencing. Can be accurate to h^2, h^4, or h^6.
    Reduces axis length by "accuracy" amount, except for zero version (special).
    See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
        * Check out that fancy recursion!
        * Uses progressively lower-accuracy methods for edges to preserve shape.
    """
    # Stuff for unequally-spaced data
    # try: x[0]
    # except TypeError:
    #     if x.ndim>1: # if want x interpreted as vector
    #         xaxis = axis
    #     else:
    #         xaxis = 0
    # else:
    #     x = np.tile(np.array([x]), y.shape[axis])
    # if x.shape[xaxis] != y.shape[axis]: # allow broadcasting rules to be used along other axes
    #     raise ValueError('x and y dimensions do not match along derivative axis.')
    # y = _permute(y, axis)
    # x = _permute(x, xaxis)
    # Simple Euler scheme
    # y = np.rollaxis(y, axis, y.ndim)
    ldiff, rdiff = (), ()
    if keepedges:
        keepleft = keepright = True
    y = np.array(y) # for safety
    y = _permute(y, axis)
    if accuracy==0:
        diff = (y[...,1:]-y[...,:-1])/h # keepleft and keepright are immaterial
    elif accuracy==2:
        diff = (1/2)*(y[...,2:]-y[...,:-2])/h
        if keepleft:
            ldiff = deriv1(h, y[...,:2], axis=-1, keepleft=True, accuracy=0), # one-tuple
        if keepright:
            rdiff = deriv1(h, y[...,-2:], axis=-1, keepright=True, accuracy=0),
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    elif accuracy==4:
        diff = (1/12)*(-y[...,4:] + 8*y[...,3:-1]
                - 8*y[...,1:-3] + y[...,:-4])/h
        if keepleft:
            ldiff = deriv1(h, y[...,:3], axis=-1, keepleft=True, accuracy=2), # one-tuple
        if keepright:
            rdiff = deriv1(h, y[...,-3:], axis=-1, keepright=True, accuracy=2),
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    elif accuracy==6:
        diff = (1/60)*(y[...,6:] - 9*y[...,5:-1] + 45*y[...,4:-2]
                - 45*y[...,2:-4] + 9*y[...,1:-5] - y[...,:-6])/h
        if keepleft:
            ldiff = deriv1(h, y[...,:5], axis=-1, keepleft=True, accuracy=4), # one-tuple
        if keepright:
            rdiff = deriv1(h, y[...,-5:], axis=-1, keepright=True, accuracy=4),
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    else:
        raise ValueError('Invalid accuracy; for now, choose form O(h^2), O(h^4), or O(h^6).')
    # return np.rollaxis(diff, y.ndim-1, axis)
    return _unpermute(diff, axis)

def deriv2(h, y, axis=0, accuracy=2, keepleft=False, keepright=False, keepedges=False):
    """
    Second order finite differencing. Can be accurate to h^2, h^4, or h^6.
    See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
        * Here, since there is no comparable midpoint-2nd derivative, need to
          just pad endpoints with the adjacent derivatives.
        * Again, check out that fancy recursion!
    Reduces axis length by "accuracy" amount, except for zero version (special).
    """
    # Simple Euler scheme
    # y = np.rollaxis(y, axis, y.ndim)
    ldiff, rdiff = (), ()
    if keepedges:
        keepleft = keepright = True
    y = np.array(y) # for safety
    y = _permute(y, axis)
    if accuracy==2:
        diff = (y[...,2:] - 2*y[...,1:-1] + y[...,:-2])/h**2
        if keepleft: # just append the leftmost 2nd deriv
            ldiff = diff[...,:1],
        if keepright: # just append the rightmost 2nd deriv
            rdiff = diff[...,-1:],
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    elif accuracy==4:
        diff = (1/12)*(-y[...,4:] + 16*y[...,3:-1]
                - 30*y[...,2:-2] + 16*y[...,1:-3] - y[...,:-4])/h**2
        if keepleft:
            ldiff = deriv2(h, y[...,:3], axis=-1, keepleft=True, accuracy=2),
        if keepright:
            rdiff = deriv2(h, y[...,-3:], axis=-1, keepright=True, accuracy=2),
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    elif accuracy==6:
        diff = (1/180)*(2*y[...,6:] - 27*y[...,5:-1] + 270*y[...,4:-2]
                - 490*y[...,3:-3] + 270*y[...,2:-4] - 27*y[...,1:-5] + 2*y[...,:-6])/h**2
        if keepleft:
            ldiff = deriv2(h, y[...,:5], axis=-1, keepleft=True, accuracy=4),
        if keepright:
            rdiff = deriv2(h, y[...,-5:], axis=-1, keepright=True, accuracy=4),
        diff = np.concatenate((*ldiff, diff, *rdiff), axis=-1)
    else:
        raise ValueError('Invalid accuracy; for now, choose form O(h^2), O(h^4), or O(h^6).')
    # return np.rollaxis(diff, y.ndim-1, axis)
    return _unpermute(diff, axis)

def deriv_uneven(x, y, axis=0, keepedges=False):
    """
    Central numerical differentiation, uneven/even spacing.
    Equation: (((x1-x0)/(x2-x1))(y2-y1) + ((x2-x1)/(x1-x0))(y1-y0)) / (x2-x0)
        * reduces to standard (y2-y0)/(x2-x0) for even spcing, and for uneven
          weights the slope closer to center point more heavily
        * want weighted average of forward/backward Euler, with weights 1 minus
          percentage of total x2-x0 interval
    Reduces axis length by 2.
    """
    # Preliminary stuff
    x, y = np.array(x), np.array(y) # precaution
    xaxis = (axis if x.ndim>1 else 0) # if want x interpreted as vector
    if x.shape[xaxis] != y.shape[axis]: # allow broadcasting rules to be used along other axes
        raise ValueError('x and y dimensions do not match along derivative axis.')
    x, y = _permute(x, xaxis), _permute(y, axis)
    # Formulation from stackoverflow, shown to be equivalent to the
    # one referenced below, but use x's instead of h's, and don't separte out terms
    # Original from this link: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
    x0, x1, x2 = x[...,:-2], x[...,1:-1], x[...,2:]
    y0, y1, y2 = y[...,:-2], y[...,1:-1], y[...,2:]
    h1, h2 = x1-x0, x2-x1 # the x-steps
    # f = (x2 - x1)/(x2 - x0)
    # diff = (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0) # version 1
    diff = -h2*y0/(h1*(h1+h2)) - (h1-h2)*y1/(h1*h2) + h1*y2/(h2*(h1+h2))
    if keepedges: # pad with simple differences on edges
        bh = np.diff(x[...,:2], axis=-1)
        eh = np.diff(x[...,-2:], axis=-1)
        bdiff = deriv1(bh, y[...,:2], axis=-1, keepedges=True, accuracy=0)
        ediff = deriv1(eh, y[...,-2:], axis=-1, keepedges=True, accuracy=0)
        diff = np.concatenate((bdiff, diff, ediff), axis=-1)
    return _unpermute(diff, axis)

def deriv1_uneven(*args, **kwargs):
    """
    Defined for name consistency.
    """
    return deriv_uneven(*args, **kwargs)

def deriv2_uneven(x, y, axis=0, keepedges=False): # alternative
    """
    Second derivative adapted from Euler's method using the same technique as above.
    """
    # Preliminary stuff
    x, y = np.array(x), np.array(y) # precaution
    xaxis = (axis if x.ndim>1 else 0) # if want x interpreted as vector
    if x.shape[xaxis] != y.shape[axis]: # allow broadcasting rules to be used along other axes
        raise ValueError('x and y dimensions do not match along derivative axis.')
    x, y = _permute(x, xaxis), _permute(y, axis)
    # Formulation from this link: https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/#comments
    # Identical to this link: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
    x0, x1, x2 = x[...,:-2], x[...,1:-1], x[...,2:]
    y0, y1, y2 = y[...,:-2], y[...,1:-1], y[...,2:]
    h1, h2, H = x1-x0, x2-x1, x2-x0 # the x-steps
    # diff = 2*((x2-x1)*y0 - (x2-x0)*y1 + (x1-x0)*y2) / ((x2-x1)*(x1-x0)*(x2-x0)) # version 1
    # diff = 2*(y0/((x1-x0)*(x2-x0)) - y1/((x2-x1)*(x1-x0)) + y2/((x2-x1)*(x2-x0))) # version 2
    diff = 2*(h2*y0 - H*y1 + h1*y2)/(h1*h2*H) # version 3
    if keepedges: # need 3 points for 2nd derivative; can only pad edges with the nearest 2nd derivs
        diff = np.concatenate((diff[...,:1], diff, diff[...,-1:]), axis=-1)
    return _unpermute(diff, axis)

def deriv3_uneven(x, y, axis=0, keepedges=False): # alternative
    """
    Second derivative adapted from Euler's method using the same technique as above.
    """
    # Preliminary stuff
    x, y = np.array(x), np.array(y) # precaution
    xaxis = (axis if x.ndim>1 else 0) # if want x interpreted as vector
    if x.shape[xaxis] != y.shape[axis]: # allow broadcasting rules to be used along other axes
        raise ValueError('x and y dimensions do not match along derivative axis.')
    x, y = _permute(x, xaxis), _permute(y, axis)
    # Formulation from the same PDF shown above
    # First 4-point formula, which is uncentered and weird, so don't use it
    # x0, x1, x2, x3 = x[...,:-3], x[...,1:-2], x[...,2:-1], x[...,3:]
    # y0, y1, y2, y3 = y[...,:-3], y[...,1:-2], y[...,2:-1], y[...,3:]
    # h1, h2, h3, H = x1-x0, x2-x1, x3-x2, x3-x0 # the x-steps
    # diff = (6/h1)*(-y0/(h1*(h1+h2)*H) + y1/(h1*h2*(h1+h3)) - y2/((h1+h2)*h2*h3) + y3/(H*(h2+h3)*h))
    # Now the 5-point formula, evaluating on the center points
    # * Changed second line from h1+h3 to h2+h3; this was just by eyeballing the
    #   coefficients and trying to enforce symmetry.
    # * Turns out the *paper is incorrect*; original version leads to weird derivative,
    #   but by changing that term the results are as expected.
    x0, x1, x2, x3, x4 = x[...,:-4], x[...,1:-3], x[...,2:-2], x[...,3:-1], x[...,4:]
    y0, y1, y2, y3, y4 = y[...,:-4], y[...,1:-3], y[...,2:-2], y[...,3:-1], y[...,4:]
    h1, h2, h3, h4 = x1-x0, x2-x1, x3-x2, x4-x3
    H1, H2, H = h1+h2+h3, h2+h3+h4, h1+h2+h3+h4 # cleaned up their notation a bit; now
        # just replace H2 with H, and replace H2-h1 with H2; this makes more sense
    # for h in h1,h2,h3,h4,H1,H2,H: print(h.min(),h.max())
    # diff = (-.5*y0 + y1 - y3 + .5*y4)/(50e2**3) # Euler method; result is actually normal
    diff = 6*((h2-2*h3-h4)*y0/(h1*(h1+h2)*H1*H) \
            - (h1+h2-2*h3-h4)*y1/(h1*h2*(h2+h3)*H2) \
            + (h1+2*h2-2*h3-h4)*y2/((h1+h2)*h2*h3*(h3+h4)) \
            - (h1+2*h2-h3-h4)*y3/(H1*(h2+h3)*h3*h4) \
            + (h1+2*h2-h3)*y4/(H*H2*(h3+h4)*h4)) # holy shitballs
    if keepedges: # need 5 points for 3rd derivative; can only pad edges with the nearest 3rd derivs
        diff = np.concatenate((diff[...,:1], diff[...,:1], diff, diff[...,-1:], diff[...,-1:]), axis=-1)
    return _unpermute(diff, axis)

def diff(x, y, axis=0):
    """
    Trivial differentiation onto half levels.
    Reduces axis length by 1.
    """
    # Preliminary stuff
    if x.ndim>1: # if want x interpreted as vector
        xaxis = axis
    else:
        xaxis = 0
    if x.shape[xaxis] != y.shape[axis]: # allow broadcasting rules to be used along other axes
        raise ValueError('x and y dimensions do not match along derivative axis.')
    # y = np.rollaxis(y, axis, y.ndim) # broadcasting rules will then help us out
    # x = np.rollaxis(x, xaxis, x.ndim)
    y = _permute(y, axis)
    x = _permute(x, xaxis)
    # Get forward difference
    return _unpermute((y[...,1:] - y[...,:-1])/(x[...,1:] - x[...,:-1]), axis)

def laplacian(lon, lat, data, accuracy=4):
    """
    Get Laplacian over geographic grid.
    Input is longitude, latitude, and data in any units.
    Equation: del^2 = (1/a^2*cos^2(phi))(d/dtheta)^2 +
                (1/a^2*cos(phi))(d/dphi)(cos(phi)*d/dphi)
        This is the geographic coordinate version.
    """
    # Setup
    npad = accuracy//2 # need +/-1 for O(h^2) approx, +/-2 for O(h^4), etc.
    data = geopad(lon, lat, data, nlon=npad, nlat=npad)[2] # pad lons/lats
    phi, theta = lat*np.pi/180, lon*np.pi/180 # from north pole
    # Execute
    h_phi, h_theta = abs(phi[2]-phi[1]), abs(theta[2]-theta[1]) # avoids weird edge cases
    phi = phi[None,...]
    for i in range(2,data.ndim):
        phi = phi[...,None]
    laplacian = ( # below avoids repeating finite differencing scheme
            (1/(const.a**2 * np.cos(phi)**2)) * deriv2(h_theta, data[:,npad:-npad,...], axis=0, accuracy=accuracy) # unpad latitudes
            + (-np.tan(phi)/(const.a**2)) * deriv1(h_phi, data[npad:-npad,...], axis=1, accuracy=accuracy) # unpad longitudes
            + (1/(const.a**2)) * deriv2(h_phi, data[npad:-npad,...], axis=1, accuracy=accuracy) # unpad longitudes
            ) # note, no axis rolling required here; the deriv schemes do that
    return laplacian

#------------------------------------------------------------------------------
# Wave activity stuff
#------------------------------------------------------------------------------
def eqlat(lon, lat, q, skip=10, sigma=None, fix=False): #n=1001, skip=10):
    """
    Get series of equivalent latitudes corresponding to PV levels evenly
    sampled from distribution of sorted PVs. Solves the equation...
        Area == integral_(-pi/2)^x a^2*2*pi*cos(phi) dphi
    ...which is the area occupied by PV zonalized below a given contour.
    Returns:
        * band, the equivalent latitude
        * qband, its associated Q threshold
        * w, grid weights (for future use)
    """
    # Initial stuff
    areas = Properties(lon, lat).areas
        # delivers grid areas, as function of latitude

    # Flatten
    q, shape = _flatten(q, end=2) # gives current q, and former shape

    # And consider mass-weighting this stuff
    if sigma is not None:
        mass = _flatten(sigma, end=2)*areas[...,None]
            # note that at least singleton dimension is added
        masscum = mass.cumsum(axis=1).sum(axis=0, keepdims=True)
            # cumulative mass from pole

    # Determing Q contour values for evaluating eqlat (no interpolation; too slow)
    K = q.shape[-1] # number of extra dims
    N = (np.prod(shape[:2])-1)//skip + 1 # e.g. // counts number of complete length-<skip> blocks after index 0, then add 0 position
    offset = (np.prod(shape[:2]) % skip)//2 # want to center the list if possible; e.g. [0,1,2,3,4,5,6,7,8], skip=3, will give [1,4,7]
    # N = N-2 # want to eliminate start/end, in case? no should be fine

    # Solve equivalent latitudes; options include mass-weighted solution with sigma, or area-weighted
    bands, q_bands = np.empty((1, N, K)), np.empty((1, N, K))
    for k in range(K): # iterate through extra dimensions
        # Iterate through Q contours
        q_bands[0,:,k] = np.sort(q[:,:,k], axis=None)[offset::skip] # test q-values
        for n in range(N):
            f = (q[:,:,k] <= q_bands[0,n,k]) # filter
            if sigma is None: # normal weighting
                # Get sine of eqlat, and correct for rounding errors
                sin = areas[f].sum()/(2*np.pi*const.a**2)-1
                if sin>1: sin = 1
                if sin<-1: sin = -1
                bands[0,n,k] = np.arcsin(sin)*180/np.pi
            else: # mass weighting
                # Interpolate to latitude of mass contour
                massk, masscumk = mass[:,:,k], masscum[:,:,k].squeeze() # latter is coordinate
                mass = massk[f].sum() # total mass below Q
                bands[0,n,k] = np.interp(mass, masscumk, lat) # simple interpolation to one point

    # Reshape, and return
    return _unflatten(bands, shape, end=2), _unflatten(q_bands, shape, end=2)

def waqlocal(lon, lat, q,
        nh=True, skip=10):
    """
    Get local wave activity measure.
    Input...
        * lon, grid longitudes
        * lat, grid latitudes
        * q, the PV stuff
        * skip, the interval of sorted q you choose (passed to eqlat)
    """
    # Grid considerations
    if nh:
        lat, q = -np.flipud(lat), -np.flip(q, axis=1)
            # negated q, so monotonically increasing "northward"
    grid = Properties(lon, lat)
    areas, dphi, phib = grid.areas, grid.dphi, grid.phib
    integral = const.a*phib[None,:]

    # Flatten (eqlat can do this, but not necessary here)
    q, shape = _flatten(q, end=2) # flatten

    # Get equivalent latiitudes
    bands, q_bands = eqlat(lon, lat, q, skip=skip) # note w is just lonbylat
    L, M, N, K = q.shape[0], q.shape[1], bands.shape[1], q.shape[-1] # number of eqlats, number of extra dims

    # Get local wave activity measure, as simple line integrals
    waq = np.empty((L, M, K))
    percent = 0
    for k in range(K):
        if (k/K)>(.01*percent):
            print('%d%% finished...' % (100*k/K,))
            percent = percent+10
        # Loop through each contour
        waq_k = np.empty((L,N))
        for n in range(N):
            # Setup, large areas
            band = bands[0,n,k]*np.pi/180
            if np.isnan(band): # happens if contours intersect at edge
                waq_k[:,n] = np.nan
            else:
                anom = q[:,:,k] - q_bands[0,n,k]
                f_pos = (anom >= 0) & (phib[None,1:] < band) # high anomalies at low lat (below top graticule)
                f_neg = (anom < 0) & (phib[None,:-1] >= band) # low anomalies at high lat (above bottom graticule)
                # See if band is sandwiched between latitudes
                mid = np.where((phib[:-1] <= band) & (phib[1:] > band))[0] # want scalar id (might be zero)
                if mid.size>0:
                    f_pos_mid = (anom[:,mid] >= 0) # longitudes where positive
                    p_int, m_int = const.a*(band-phib[mid]), const.a*(phib[mid+1]-band)
                        # partial integrals, positive and negative
                for l in range(L):
                    # Get individual integral
                    integral_pos = (anom[l,f_pos[l,:]]*integral[:,f_pos[l,:]]).sum()
                    integral_neg = -(anom[l,f_neg[l,:]]*integral[:,f_neg[l,:]]).sum() # minus a negative
                    if mid.size>0:
                        if f_pos_mid[l]: # if positive at this latitude, we add anomaly
                            integral_extra = anom[l,mid]*p_int
                        else: # else, subtract it
                            integral_extra = -anom[l,mid]*m_int
                    else:
                        integral_extra = 0
                    # Put it all together
                    waq_k[l,n] = integral_pos + integral_neg + integral_extra # no normalization here
        # Interpolate
        for l in range(L):
            waq[l,:,k] = np.interp(lat, bands[0,:,k], waq_k[l,:])

    # Return
    if nh: waq = np.flip(waq, axis=1)
    return _unflatten(waq, shape)

def waq(lon, lat, q, sigma=None, omega=None,
        nh=True, skip=10): #, ignore=None): #N=1001, ignore=None):
    """
    Get finite-amplitude wave activity.
    Input...
        * lon, grid longitudes
        * lat, grid latitudes
        * q, the PV quantity
        * sigma (optional), the instability
        * omega (optional), the quantity being integrated; note the isentropic mass equation
        * skip, the interval of sorted q you choose (passed to eqlat)
        * nh (bool)
    """
    # Grid considerations
    if nh:
        lat, q = -np.flipud(lat), -np.flip(q, axis=1)
        if omega is not None: omega = -np.flip(omega, axis=1)
        if sigma is not None: sigma = np.flipd(sigma, axis=1)
        # negated q/omega, so monotonically increasing "northward"
    grid = Properties(lon, lat)
    areas, dphi, phib = grid.areas, grid.dphi, grid.phib

    # Flatten (eqlat can do this, but not necessary here)
    q, shape = _flatten(q, end=2) # flatten
    if omega is not None: omega, _ = _flatten(omega, end=2) # flatten
    if sigma is not None: sigma, _ = _flatten(sigma, end=2)

    # Get equivalent latiitudes
    bands, q_bands = eqlat(lon, lat, q, sigma=sigma, skip=skip) # note w is just lonbylat
        # will infer area weights, to get equivalent latitude
    M, N, K = q.shape[1], bands.shape[1], q.shape[-1] # number of eqlats, number of extra dims

    # Get activity
    waq = np.empty((1, M, K))
    percent = 0
    for k in range(K):
        if (k/K)>(.01*percent):
            print('%d%% finished...' % (percent,))
            percent = percent+10
        # Loop through each contour
        waq_k = np.empty(N)
        for n in range(N): #i, bandki in enumerate(bandk[0,:,k]): #, q_bandki) in enumerate(zip(bandk, q_bandk)):
            # First, main blocks
            band = bands[0,n,k]*np.pi/180
            if np.isnan(band):
                waq_k[n] = np.nan
            else:
                # anom = q[:,:,k] - q_bands[0,n,k]
                qk, Qk = q[:,:,k], q_bands[0,n,k] # should be identical, since positive/negative regions have same area by construction
                if omega is None:
                    qint = q[:,:,k] # the thing being integrated
                else:
                    qint = omega[:,:,k]
                # f_pos = (anom >= 0) & (phib[None,1:] < band) # high anomalies at low lat (below top graticule)
                # f_neg = (anom < 0) & (phib[None,:-1] >= band) # low anomalies at high lat (above bottom graticule)
                # integral = (anom[f_pos]*areas[f_pos]).sum() - (anom[f_neg]*areas[f_neg]).sum() # minus a negative
                f_pos = (qk >= Qk) & (phib[None,1:] < band) # high anomalies at low lat (below top graticule)
                f_neg = (qk < Qk) & (phib[None,:-1] >= band) # low anomalies at high lat (above bottom graticule)
                integral = (qint[f_pos]*areas[f_pos]).sum() - (qint[f_neg]*areas[f_neg]).sum() # minus a negative
                # Next, account for tiny pieces along equivalent latitude cells
                mid = np.where((phib[:-1] <= band) & (phib[1:] > band))[0] # want scalar id
                try: mid = mid[0]
                except IndexError:
                    integral_extra = 0
                else:
                    # f_pos_mid = (anom[:,mid] >= 0) # longitudes where positive
                    f_pos_mid = (qk[:,mid] >= Qk) # longitudes where positive
                    p_dphi, m_dphi = (
                            np.cos((band+phib[mid])/2)*(band-phib[mid]), # positive, low lat
                            np.cos((band+phib[mid+1])/2)*(phib[mid+1]-band) # negative, high lat
                            )
                    integral_extra = (
                            qint[f_pos_mid,mid].sum()*(areas[mid]*m_dphi/dphi[mid])
                            - qint[~f_pos_mid,mid].sum()*(areas[mid]*p_dphi/dphi[mid])
                            )
                # Put it all together
                waq_k[n] = (integral + integral_extra)/(2*np.pi*const.a*np.cos(band))
        # Interpolate
        nanfilt = np.isnan(waq_k)
        if sum(~nanfilt)==0:
            print('Warning: no valid waqs calculated for k %d.' % k)
            waq[0,:,k] = np.nan
        else:
            waq[0,:,k] = np.interp(grid.latc, bands[0,~nanfilt,k], waq_k[~nanfilt])

    # Return
    if nh: waq = np.flip(waq, axis=1)
    return _unflatten(waq, shape)

#------------------------------------------------------------------------------
# Changing phase space (to EOFs, spectral decomposition, etc.)
#------------------------------------------------------------------------------
def eof(data, neof=5):
    """
    Calculates the temporal EOFs, using most efficient method.
    """
    return

def autocorr():
    """
    Gets the autocorrelation spectrum at successive lags.
    """
    return

def rednoisefit():
    """
    Returns a best-fit rednoise autocorrelation spectrum.
    """
    return

def lowpass(x, k=4, axis=-1): #n=np.inf, kmin=0, kmax=np.inf): #, kscale=1, krange=None, k=None):
    """
    Extracts the time series associated with cycle, given by the first k
    Fourier harmonics for the time series.
    Does not apply any windowing.
    """
    # Naively remove certain frequencies
    # p = np.abs(fft)**2
    # f = np.where((freq[1:]>=kmin) | (freq[1:]<=kmax))
    #     # should ignore first coefficient, the mean
    # if n==np.inf: fremove = f
    # else: fremove = f[np.argpartition(p, -n)[-n:]]
    #         # gets indices of n largest values

    # # And filter back
    # fft[fremove+1], fft[-fremove-1] = 0+0j, 0+0j
    #     # fft should be symmetric, so remove locations at corresponding negative freqs

    # Get fourier transform
    # x = np.rollaxis(x, axis, x.ndim)
    x = _permute(x, axis)
    fft = np.fft.fft(x, axis=-1)
    # freq = np.fft.fftfreq(x.size)*scale

    # Remove the edge case frequencies
    fft[...,0] = 0
    fft[...,k+1:-k] = 0
    # return np.rollaxis(np.fft.ifft(fft).real, x.ndim-1, axis)
    return _unpermute(np.fft.ifft(fft).real, axis)
        # FFT will have some error, and give non-zero imaginary components;
        # just naively cast to real

def lanczos(alpha, J):
    """
    Lanczos filtering of data; gives an abrupt high-frequency (low wavenumber)
    cutoff at omega = alpha*pi, the number of datapoints needed.
    """
    C0 = alpha # integral of cutoff-response function is alpha*pi/pi
    Ck = np.sin(alpha*np.pi*np.arange(1,J+1))*(1/(np.pi*np.arange(1,J+1)))
    Cktilde = Ck*np.sin(np.pi*np.arange(1,J+1)/J)/(np.pi*np.arange(1,J+1)/J)
    filt = np.concatenate((np.flipud(Cktilde), np.array([C0]), Cktilde))
    return filt/filt.sum()
    # for j,J in enumerate(Jsamp):
    #     pass
    # R = lambda Cfunc, omega: C0 + np.sum(
    #         Cfunc(alpha)*np.cos(omega*np.arange(1,J+1))
    #         ) # C_tau * cos(omega*tau)
    # omega = np.linspace(0,np.pi,1000)
    # Romega = np.empty(omega.size)
    # Romegatilde = np.empty(omega.size)
    # for i,o in enumerate(omega): Romega[i] = R(Ck, o)
    # for i,o in enumerate(omega): Romegatilde[i] = R(Cktilde, o)
    # a1.plot(omega, Romega, color=colors[j], label=('J=%d' % J))
    # a2.plot(omega, Romegatilde, color=colors[j], label=('J=%d' % J))
    # return

def butterworth(*args, **kwargs):
    """
    Return Butterworth filter.
    Wraps around the builtin scipy method.
    """
    b, a = signal.butter(*args, **kwargs)
    return b/a # gives numerate/demoninator of coefficients; we just want floating points

def window(wintype, M=100):
    """
    Retrieves weighting function window.
    """
    # Prepare window
    if wintype=='boxcar':
        win = np.ones(M)
    elif wintype=='hanning':
        win = np.hanning(M) # window
    elif wintype=='hamming':
        win = np.hamming(M)
    elif wintype=='blakman':
        win = np.blackman(M)
    elif wintype=='kaiser':
        win = np.kaiser(M, param)
    elif wintype=='lanczos':
        win = lanczos(M, param)
    elif wintype=='butterworth':
        win = butterworth(M, param)
    else:
        raise ValueError('Unknown window type: %s' % (wintype,))
    return win/win.sum()

def spectrum(x, M=72, wintype='boxcar', param=None, axis=-1):
    """
    Gets the spectral decomposition for particular windowing technique.
    """
    # Initital stuff; get integer number of half-overlapped windows
    N = x.size
    pm = M//2
    Nround = pm*(N//pm)
    x = x[:Nround]
    if N-Nround>0: print(f'Points removed: {N-Nround:d}.')

    # Get copsectrum, quadrature spectrum, and powers for each window
    win = window(wintype, M)
    loc = np.arange(pm, Nround-pm+pm//2, pm) # jump by half window length
    Cx = np.empty((loc.size, pm)) # have half/window size number of freqs
    for i,l in enumerate(loc):
        Cx[i,:] = np.abs(np.fft.fft(win*signal.detrend(x[l-pm:l+pm]))[:pm])**2
        # numpy fft gives power A+Bi, so want sqrt(A^2 + B^2)/2 for +/- wavenumbers
    freq = np.fft.fftfreq(M)[:pm] # frequency
    return freq, Cx.mean(axis=0)

def cspectrum(x, y, M=72, wintype='boxcar', param=None, centerphase=np.pi):
    """
    Calculates the cross spectrum for particular windowing technique; kwargs
    are passed to scipy.signal.cpd
    """
    # Iniital stuff; get integer number of half-overlapped windows
    N = x.size
    pm = M//2
    Nround = pm*(N//pm)
    x, y = x[:Nround], y[:Nround]
    if N-Nround>0: print(f'Points removed: {N-Nround:d}.')

    # Get copsectrum, quadrature spectrum, and powers for each window
    win = window(wintype, M)
    loc = np.arange(pm, Nround-pm+pm//2, pm) # jump by half window length
    shape = (loc.size, pm) # have half window size number of freqs
    Fxx, Fyy, CO, Q = np.empty(shape), np.empty(shape), np.empty(shape), np.empty(shape)
    for i,l in enumerate(loc):
        Cx = np.fft.fft(win*signal.detrend(x[l-pm:l+pm]))[:pm]
        Cy = np.fft.fft(win*signal.detrend(y[l-pm:l+pm]))[:pm]
        Fxx[i,:] = np.abs(Cx)**2
        Fyy[i,:] = np.abs(Cy)**2
        CO[i,:] = Cx.real*Cy.real + Cx.imag*Cy.imag
        Q[i,:] = Cx.real*Cy.imag - Cy.real*Cx.imag

    # Get average cospectrum and other stuff, return
    Fxx, Fyy, CO, Q = Fxx.mean(0), Fyy.mean(0), CO.mean(0), Q.mean(0)
    Coh = (CO**2 + Q**2)/(Fxx*Fyy) # coherence
    p = np.arctan2(Q, CO) # phase
    p[p >= centerphase+np.pi] -= 2*np.pi
    p[p < centerphase-np.pi] += 2*np.pi
    freq = np.fft.fftfreq(M)[:pm] # frequency
    return freq, Coh, p

def autospectrum():
    """
    Uses scipy.signal.welch windowing method to generate an estimate of the spectrum.
    """
    return

def autocspectrum():
    """
    Uses scipy.signal.cpd automated method to generate cross-spectrum estimate.
    """
    return

#-------------------------------------------------------------------------------
# TODO: Finish harvesting this original function written for the Objective
# Analysis assignment. Consider deleting it.
#-------------------------------------------------------------------------------
def spectral(fnm, nm, data, norm=True, win=501,
            freq_scale=1, scale='days',
            xlog=True, ylog=False, mc='k',
            xticks=None,
            rcolors=('C3','C6'), pcolors=('C0','C1'), alpha=0.99, marker=None,
            xlim=None, ylim=None, # optional override
            linewidth=1.5,
            red_discrete=True, red_contin=True, manual=True, welch=True):
    '''
    Spectral transform function. Needs work; was copied from OA
    assignment and right now just plots a bunch of stuff.
    '''
    # Iniital stuff
    N = len(data)
    pm = int((win-1)/2)
    fig, a = plt.subplots(figsize=(15,5))

    # Confidence intervals
    dof_num = 1.2*2*2*(N/(win/2))
    trans = 0.5
    exp = 1
    F99 = stats.f.ppf(1-(1-alpha)**exp,dof_num,1000)
    F01 = stats.f.ppf((1-alpha)**exp,dof_num,1000)
    print('F stats:',F01,F99)
    rho = np.corrcoef(data[1:],data[:-1])[0,1]
    kr = np.arange(0,win//2+1)
    fr = freq_scale*kr/win

    def xtrim(f):
        if xlim is None:
            return np.ones(f.size, dtype=bool)
        else:
            return ((f>=xlim[0]) & (f<=xlim[-1]))

    # Power spectra
    if manual:
        label = 'power spectrum'
        if welch: label = 'manual method'
        # Now, manual method with proper overlapping etc.
        if False:
            data = data[:(N//pm)*pm]
        loc = np.linspace(pm,N-pm-1,2*int(np.round(N/win))).round().astype(int) # sample loctaions
        han = np.hanning(win)
        han = han/han.sum()
        phi = np.empty((len(loc),win//2))
        for i,l in enumerate(loc):
            pm = int((win-1)/2)
            C = np.fft.fft(han*signal.detrend(data[l-pm:l+pm+1]))
            phii = np.abs(C)**2/2
            phii = 2*phii[1:win//2+1]
            phi[i,:] = phii
        phi = phi.mean(axis=0)
        print('phi sum:',phi.sum())
        f = np.fft.fftfreq(win)[1:win//2+1]*freq_scale
        if norm: phi = phi/phi.sum()
        f, phi = f[xtrim(f)], phi[xtrim(f)] # trim
        a.plot(f, phi, label=label,
               mec=mc, mfc=mc, mew=linewidth,
               marker=marker, color=pcolors[0], linewidth=linewidth)
        if xlim is None: xlim = ((f*freq_scale).min(), (f*freq_scale).max())
        if ylim is None: ylim = ((phi.min()*0.95, phi.max()*1.05))

    if welch:
        label = 'power spectrum'
        if manual: label = 'welch method'
        # Welch
        fw, phi_w = signal.welch(data, nperseg=win, detrend='linear', window='hanning', scaling='spectrum',
                              return_onesided=False)
        fw, phi_w = fw[1:win//2+1]*freq_scale, phi_w[1:win//2+1]
        if norm: phi_w = phi_w/phi_w.sum()
        fw, phi_w = fw[xtrim(fw)], phi_w[xtrim(fw)] # trim
        print('phiw sum:',phi_w.sum())
        a.plot(fw, phi_w, label=label,
              mec=mc, mfc=mc, mew=linewidth,
               marker=marker, color=pcolors[-1], linewidth=linewidth)
        if xlim is None: xlim = ((fw).min(), (fw).max())
        if ylim is None: ylim = (phi_w.min()*0.95, phi_w.max()*1.05)

    # Best fit red noise spectrum
    if red_discrete:
        print('Autocorrelation',rho)
        phi_r1 = (1-rho**2)/(1+rho**2-2*rho*np.cos(kr*np.pi/(win//2)))
        print('phi_r1 sum:',phi_r1.sum())
        if norm: phi_r1 = phi_r1/phi_r1.sum()
        frp, phi_r1 = fr[xtrim(fr)], phi_r1[xtrim(fr)]
        a.plot(fr[xtrim(fr)], phi_r1, label=r'red noise, $\rho(\Delta t)$',
               marker=None, color=rcolors[0], linewidth=linewidth)
        a.plot(frp, phi_r1*F99, linestyle='--',
               marker=None, alpha=trans, color=rcolors[0], linewidth=linewidth)

    # Alternate best fit
    if red_contin:
        Te = -1/np.log(rho)
        omega = (kr/win)*np.pi*2
        phi_r2 = 2*Te/(1+(Te**2)*(omega**2))
        print('phi_r2 sum:',phi_r2.sum())
        if norm: phi_r2 = phi_r2/phi_r2.sum()
        frp, phi_r2 = fr[xtrim(fr)], phi_r2[xtrim(fr)]
        a.plot(frp, phi_r2, label=r'red noise, $T_e$',
               marker=None, color=rcolors[1], linewidth=linewidth)
        a.plot(frp, phi_r2*F99, linestyle='--',
               marker=None, alpha=trans, color=rcolors[-1], linewidth=linewidth)
    # Variance
    print('true variance:',data.std()**2)
    # Figure formatting
    a.legend()
    if ylog:
        a.set_yscale('log')
    if xlog:
        a.set_xscale('log')
    a.set_title('%s power spectrum' % nm)
    a.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%5.3g'))
    a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    if xticks is None: xticks = a.get_xticks()
    my.format(a, xlabel=('frequency (%s${}^{-1}$)' % scale), ylabel='proportion variance explained',
             xlim=xlim, ylim=ylim, xticks=xticks)
    suffix = 'pdf'
    fig.savefig('a5_' + fnm + '.' + suffix, format=suffix, dpi='figure')
    plt.show()

#------------------------------------------------------------------------------
# TODO: Physical meteorological quantities
# This sections is *very limited*, and should consider deleting it
# Should probably be replaced with *climate.py* stuff
#------------------------------------------------------------------------------
def pt(lev, T, axis=2):
    print(lev.shape, T.shape)
    """
    Get potential temperature from T; by default, assume the height dimension
    is axis=2. Axis is for the instance where we have non-vector T input.
    Each level should be in mb/hpa.
    """
    if T.ndim!=1: # vector
        axis = axis % T.ndim # modify if, for example, axis=-1
        for i in range(axis):
            lev = lev[None,...]
        for i in range(T.ndim-1-axis):
            lev = lev[...,None]
    return T*(const.p0/lev)**const.kappa # ...each p here should be in mb/hPa

def absvo(lat, Zeta):
    """
    Gets absolute vorticity from latitude and vorticity.
    """
    # Get Coriolis force
    f = 2*const.Omega*np.sin(lat*np.pi/180)[None,:,None]
    for i in range(3, Zeta.ndim): f = f[...,None] # add extra dimensions

    # Final answer; account for levels lost in finite difference
    return (Zeta + f)

def pv(lat, theta, P, Zeta, uneven=True):
    """
    Gets Ertel's PV on theta-surface, from the pressure level info and
    vorticity... note theta should be in K.
    Returns...
        * absolute vorticity, the numerator
        * mass factor sigma, the denominator
    To get PV yourself, just divide them, but might need components.
    Already provided by ERA-Interim.
    """
    # Get dp/dtheta
    if uneven:
        sigma = -(1/const.g)*deriv_uneven(theta, P, axis=2) # use our differentiation method; reduces size along axis=2 by 2
    else:
        sigma = -(1/const.g)*diff(theta, P, axis=2)

    # Get Coriolis force
    f = 2*const.Omega*np.sin(lat*np.pi/180)[None,:,None]
    for i in range(3, P.ndim): f = f[...,None] # add extra dimensions

    # Final answer; account for levels lost in finite difference
    s = slice(1,-1) if uneven else slice(1,None)
    return (Zeta + f)[:,:,s,...], sigma

def qgpv(lon, lat, lev, TT, Phi,
        uneven=True, forward=True, fill=True, latmin=5, latmax=85):
    """
    Get QGPV in pseudo-height coordinates (cf. Nakamura and Solomon, 2010: Eq (2))
        * Input T is in K, geopotential Phi in m2/s2, pressure in mb/hPa.
        * Input uneven (bool) says if we use uneven vs. even finite differencing in height.
        * Input forward (bool) says whether we use forward or backward Euler, if not doing uneven method
        * Input fill (bool) says whether we fill points near poles/equator with NaN.
        If the former, input geopotential heights in m;
        If the latter, input geopotential (NOT heights, actual geopotential)
    Equation: qg = f + zeta + (f/rho0)(d/dz)(rho0*(theta-thetabar)/(d/dz)thetabar)
    REDUCES HEIGHT DIMENSION LENGTH BY 4, PRESERVES LON/LAT DIMENSIONS BY WRAPPING; RETURNS
    ADJUSTED METADATA.

    * Note that, with Z = -Hln(p/p0), have dZ = -H(p0/p)(dp/p0) = -Hd(ln(p)) simply
    log-derivative in height, so the (d/dz)'s in the right-hand term cancel.
    """
    # Dimensions, and preparations
    Z = -const.H*np.log(lev/const.p0)
    if not all(Tsh==psh for Tsh,psh in zip(TT.shape, Phi.shape)):
        raise ValueError('Temperature and vorticity parameter have different shapes.')
    if lev.size<=3:
        raise ValueError('Need at least three levels to take vertical derivative.')

    # The simple 1-D params
    f = 2*const.Omega*np.sin(lat*np.pi/180)[None,:,None] # goes into full qgpv formula
    rho0 = np.exp(-Z/const.H)[None,None,:] # actually is propto, but constants next to exp cancel out below
    for i in range(3, TT.ndim): # add in extra dimensions
        f, rho0 = f[...,None], rho0[...,None]

    # Calculate theta
    theta = pt(lev, TT) # simple as that

    # Start clock
    t = T.clock()

    # The theta-params
    print('Calculating global mean potential temperatures.')
    thetabar = geomean(lon, lat, theta, keepdims=True) # keep dims for broadcasting later
    t, tb = T.clock(), t # reassign t to tb "t before", get new time
    print('Global mean theta: %s seconds' % (t-tb))
    f_deriv = deriv_uneven if uneven else diff
    dthetabar_dz = f_deriv(Z, thetabar, axis=2) # use our differentiation method; reduces size along axis=2 by 2
    t, tb = T.clock(), t
    print('Vertical gradient global mean theta: %s seconds' % (t-tb))

    # The giant "stretching" term
    # First, the slice
    if uneven:
        slice_outer, slice_inner = slice(2,-2), slice(1,-1) # used deriv_uneven above
    elif forward:
        slice_outer, slice_inner = slice(None,-2), slice(None,-1) # used diff above above
    else:
        slice_outer, slice_inner = slice(2,None), slice(1,None)
    # And calculate
    h = (f/rho0)[:,:,slice_outer,...] * deriv_uneven(Z[slice_inner], (rho0*(theta-thetabar))[:,:,slice_inner,...]/dthetabar_dz, axis=2)
        # try simple differentiation instead
    t, tb = T.clock(), t
    print('Stretching term: %s seconds' % (t-tb))

    # Geostrophic relative vorticity
    # From deriv_uneven above
    if uneven:
        zetag = (1/f)*laplacian(lon, lat, Phi[:,:,2:-2,...], accuracy=2)
    # Diff above
    elif forward:
        zetag = (1/f)*laplacian(lon, lat, Phi[:,:,:-2,...], accuracy=2)
    else:
        zetag = (1/f)*laplacian(lon, lat, Phi[:,:,2:,...], accuracy=2)
    t, tb = T.clock(), t
    print('Geostrophic vorticity: %s seconds' % (t-tb))

    # Adjust for the middle ones, and return
    # if nanfill: q[:,i,...] = np.nan
    # else: q[:,i,...] = q[:,i,...].mean(axis=0, keepdims=True) # will broadcast
    q = f + zetag + h
    # q = h
    # q = f + zetag
    if fill:
        for i,l in enumerate(lat): # ...makes sense to use central latitude
            if (np.abs(l)<=latmin) or (np.abs(l)>=latmax): q[:,i,...] = q[:,i,...].mean(axis=0,keepdims=True)
    # Return, finally
    return q

