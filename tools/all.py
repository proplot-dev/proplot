#!/usr/bin/env python3
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
import scipy.optimize as optimize
from .. import const
# import time as T # have some big algorithms here; want to guage how long they take
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
def _trail_flatten(data, nflat):
    """
    Flatten *trailing* dimensions onto single dimension (useful).
    Call this to manipulate data on first dimensions, then restore with unflatten.
    * Note that numpy prod of an empty iterable will be 1; so, adds singleton dim.
    * Note we use tuple expansion of the [shape] tuple
    * Note default for numpy is is row-major.
    """
    return np.reshape(data, (*data.shape[:-nflat], np.prod(data.shape[-nflat:]).astype(int)), order='F'), data.shape

def _trail_unflatten(data, shape, nflat):
    """
    Undo action of flatten.
    Shape can be the original shape, or a new shape.
    """
    if data.shape[-1] != np.prod(shape[-nflat:]):
        raise ValueError(f'Number of trailing elements {data.shape[-1]} does not match trailing shape {shape[nflat:]:s}.')
    if not all(s1==s2 for s1,s2 in zip(data.shape[:-1], shape[:-nflat])):
        raise ValueError(f'Leading dimensions on data, {data.shape[:-1]}, do not match leading dimensions on new shape, {shape[:-nflat]}.')
    return np.reshape(data, shape, order='F')

def _lead_flatten(data, nflat):
    """
    Flatten *leading* dimensions onto single dimension.
    """
    return np.reshape(data, (np.prod(data.shape[:nflat]).astype(int), *data.shape[nflat:]), order='C') # make column major

def _lead_unflatten(data, shape, nflat):
    """
    Undo action of leadflatten.
    """
    if data.shape[0] != np.prod(shape[:nflat]):
        raise ValueError(f'Number of leading elements {data.shape[0]} does not match leading shape {shape[end:]:s}.')
    if not all(s1==s2 for s1,s2 in zip(data.shape[1:], shape[nflat:])):
        raise ValueError(f'Trailing dimensions on data, {data.shape[1:]}, do not match tailing dimensions on new shape, {shape[nflat:]}.')
    return np.reshape(data, shape, order='F')

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
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
    # Load dataset
    if not os.path.exists(filename):
        raise ValueError(f'{filename} does not exist!')
    with xr.open_dataset(filename, engine='netcdf4', cache=False, **kwargs) as file:
        #----------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
    # Fix precision of time units... some notes:
    # 1) sometimes get weird useless round-off error; convert to days, then restore to numpy datetime64[ns]
    #   because xarray seems to require it
    # 2) ran into mysterious problem where dataset could be loaded, but then
    #   COULD NOT BE SAVED because one of the datetimes wasn't serializable... this was
    #   in normal data, the CCSM4 CMIP5 results, made no sense; range was 1850-2006
    if 'time' in data.indexes:
        data['time'] = data.time.values.astype('datetime64[D]').astype('datetime64[ns]')
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
    # Make latitudes monotonic (note extracting values way way faster)
    try: data.lat.values[1]
    except IndexError:
        pass
    else:
        if data.lat.values[0]>data.lat.values[1]:
            data = data.isel(lat=slice(None,None,-1))
    #--------------------------------------------------------------------------#
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
        hours=(0,6,12,18), hour=None,
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
    #--------------------------------------------------------------------------#
    # Data stream
    stream = { # oper is original, moda is monthly mean of daily means
            'synoptic':  'oper',
            'monthly':   'moda'
            }.get(stream)
    if stream is None:
        raise ValueError('Must choose from "oper" for synoptic fields, "moda" for monthly means of daily means.')
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
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
            dates = '/'.join('%04d%02d00' % (y0 + (m0+n-1)//12, (m0+n-1)%12 + 1) for n in range(N))
        else:
            dates = '/to/'.join(d.strftime('%Y%m%d') for d in daterange) # MARS will get calendar days in range
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
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
    #--------------------------------------------------------------------------#
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
    # The astype(int) is actually super fast (ns), and below is in general
    # faster than list comprehension with container of datetime objects
    # UNIX time starts at 1970-01-01 00:00:00
    return dt.astype('datetime64[Y]').astype(np.int32)+1970

def month(dt):
    """
    Gets month from numpy datetime object (used e.g. by xarray, pandas).
    """
    # Below will convert datetime64 units from [ns] (default) to months, then spit out months relative to year
    # UNIX time starts at 1970-01-01 00:00:00
    return dt.astype('datetime64[M]').astype(np.int32)%12 + 1

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
def intersection(x, segment1, segment2, xlog=False):
    """
    Find the (first) intersection point for two line segments.
    Optionally do this in log-space for the x-axis.
    """
    #--------------------------------------------------------------------------#
    # Initial stuff
    segment1, segment2 = np.array(segment1), np.array(segment2)
    if xlog:
        func = lambda x: np.log10(x)
        ifunc = lambda x: 10**x
    else:
        func = lambda x: x
        ifunc = lambda x: x
    #--------------------------------------------------------------------------#
    # Get intersection
    diff = segment1 - segment2
    if (diff>0).all() or (diff<0).all():
        print("Warning: No intersections found.")
        return np.nan, np.nan
    idx = np.where(diff>0)[0][0]
    x, y = diff[idx-1:idx+1], func(x[idx-1:idx+1]) # two-element vectors
    px = ifunc(y[0] + (0-x[0])*((y[1]-y[0])/(x[1]-x[0])))
    x, y = y, segment2[idx-1:idx+1] # once again for the y-position; can also use segment1
    py = (y[0] + (func(px)-x[0])*((y[1]-y[0])/(x[1]-x[0])))
    return px, py

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
    if axis<0: axis = data.ndim+axis # e.g. if 3 dims, and want to axis dim -1, this is dim number 2
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
    # Coefficients are returned in reverse; e.g. 2-deg fit gives c[0]*x^2 + c[1]*x + c[2]
    slopes = coeff[0,:]
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
# Changing phase space (to EOFs, spectral decomposition, etc.)
#------------------------------------------------------------------------------
def rednoise(ntime, nsamples, a, mean=0, stdev=1, nested=False):
    """
    Creates artificial red noise time series, i.e. a weighted sum of random perturbations.
    Equation is: x(t) = a*x(t-dt) + b*eps(t)
    where a is the lag-1 autocorrelation and b is a scaling term.
     * Output will have shape ntime by nsamples.
     * Enforce that the first timestep always equals the 'starting' position.
     * Use 'nested' flag to control which algorithm to use.
    """
    #--------------------------------------------------------------------------#
    # Initial stuff
    ntime -= 1 # exclude the initial timestep
    output = np.empty((ntime+1,nsamples))
    output[0,:] = 0 # initiation
    b = (1-a**2)**0.5 # from OA class
    #--------------------------------------------------------------------------#
    # Nested loop
    for i in range(nsamples):
        eps = np.random.normal(loc=0, scale=1, size=ntime)
        for t in range(1,ntime+1):
            output[t,i] = a*output[t-1,i] + b*eps[t-1]
    #--------------------------------------------------------------------------#
    # This formula was nonsense; only works for autocorrelation 1
    # return np.concatenate((np.ones((1,nsamples))*mean,
    #   mean + np.random.normal(loc=0, scale=stdev, size=(ntime-1,nsamples)).cumsum(axis=0)), axis=0).squeeze()
    #--------------------------------------------------------------------------#
    # Trying to be fancy, just turned out fucking way slower
    # aseries = b*np.array([a**(ntime-i) for i in range(1,ntime+1)])
    # for i in range(nsamples):
    #     eps = np.random.normal(loc=0, scale=1, size=ntime)
    #     vals = eps[:,None]@aseries[None,:] # matrix for doing math on
    #     output[1:,i] = [np.trace(vals,ntime-i) for i in range(1,ntime+1)]
    return mean + stdev*output.squeeze() # rescale to have specified stdeviation/mean

def eof(data, neof=5):
    """
    Calculates the temporal EOFs, using most efficient method.
    """
    return

def autocorr(data, nlag=None, lag=None, verbose=False, axis=0):
    """
    Gets the autocorrelation spectrum at successive lags.
      * Estimator is: ((n-k)*sigma^2)^-1 * sum_i^(n-k)[X_t - mu][X_t+k - mu]
        See: https://en.wikipedia.org/wiki/Autocorrelation#Estimation
      * By default, includes lag-zero values; if user just wants single lag
        will, however, throw out those values.
    """
    #--------------------------------------------------------------------------#
    # Preparation, and stdev/means
    data = np.array(data)
    naxis = data.shape[axis] # length
    if (nlag is None and lag is None) or (nlag is not None and lag is not None):
        raise ValueError(f"Must specify *either* a lag (\"lag=x\") or range of lags (\"nlag=y\").")
    if nlag is not None and nlag>=naxis/2:
        raise ValueError(f"Lag {nlag} must be greater than axis length {naxis}.")
    if verbose:
        if nlag is None:
            print(f"Calculating lag-{lag} autocorrelation.")
        else:
            print(f"Calculating autocorrelation spectrum up to lag {nlag} for axis length {naxis}.")
    data = _permute(data, axis)
    mean = data.mean(axis=-1, keepdims=True) # keepdims for broadcasting in (data minus mean)
    var = data.var(axis=-1, keepdims=False) # this is divided by the summation term, so should have annihilated axis
    #--------------------------------------------------------------------------#
    # Loop through lags
    # Include an extra lag under some circumstances
    if nlag is None and lag==0:
        autocorrs = np.ones((*data.shape[:-1],1))
    elif nlag is None:
        autocorrs = np.sum((data[...,:-lag]-mean)*(data[...,lag:]-mean),axis=-1)/((naxis-lag)*var)
        autocorrs = autocorrs[...,None] # add dimension back in
    else:
        autocorrs = np.empty((*data.shape[:-1], nlag+1)) # will include the zero-lag autocorrelation
        autocorrs[...,0] = 1 # lag-0 autocorrelation
        for i,lag in enumerate(range(1,nlag+1)):
            autocorrs[...,i+1] = np.sum((data[...,:-lag]-mean)*(data[...,lag:]-mean),axis=-1)/((naxis-lag)*var)
    return _unpermute(autocorrs, axis)

def rednoisefit(data, nlag=None, axis=-1, lag1=False, series=False, verbose=False):
    """
    Returns a best-fit red noise autocorrelation spectrum.
    * If 'series' is True, return the red noise spectrum full series.
    * If 'series' is False, return just the *timescale* associated with the red noise spectrum.
    * If 'lag1' is True, just return the lag-1 version.
    Go up to nlag-timestep autocorrelation.
    """
    #--------------------------------------------------------------------------#
    # Initial stuff
    if nlag is None:
        raise ValueError(f"Must declare \"nlag\" argument; number of points to use for fit.")
    nflat = data.ndim-1
    data = _lead_flatten(_permute(data, axis), nflat)
    # if len(time)!=data.shape[-1]:
    #     raise ValueError(f"Got {len(time)} time values, but {data.shape[-1]} timesteps for data.")
    #--------------------------------------------------------------------------#
    # First get the autocorrelation spectrum, and flatten leading dimensions
    # Dimensions will need to be flat because we gotta loop through each 'time series' and get curve fits.
    # time = time[:nlag+1] # for later
    autocorrs = autocorr(data, nlag, axis=-1, verbose=verbose)
    #--------------------------------------------------------------------------#
    # Next iterate over the flattened dimensions, and perform successive curve fits
    ndim = data.shape[-1] if series else 1
    shape = (*data.shape[:-1], ndim) # when we perform *unflattening*
    output = np.empty((autocorrs.shape[0], ndim))
    time = np.arange(autocorrs.shape[-1]) # time series for curve fit
    for i in range(autocorrs.shape[0]): # iterate along first dimension; each row is an autocorrelation spectrum
        if lag1:
            # dt = time[1]-time[0] # use this for time diff
            p = [-1/np.log(autocorrs[i,1])] # -ve inverse natural log of lag-1 autocorrelation
        else:
            p, _ = optimize.curve_fit(lambda t,tau: np.exp(-t/tau), time, autocorrs[i,:])
        if series:
            output[i,:] = np.exp(-np.arange(ndim)/p[0]) # return the best-fit red noise spectrum
        else:
            output[i,0] = p[0] # add just the timescale
    return _unpermute(_lead_unflatten(output, shape, nflat), axis)

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

def butterworth(order, *args, **kwargs):
    """
    Return Butterworth filter.
    Wraps around the builtin scipy method.
    """
    kwargs.update({'analog':True}) # must be analog or something
    b, a = signal.butter(order-1, *args, **kwargs)
    return b/a # gives numerate/demoninator of coefficients; we just want floating points
        # an order of N will return N+1 denominator values; we fix this

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

#------------------------------------------------------------------------------#
# TODO Finish parsing this; was copied from sst project
# Gets the significance accounting for autocorrelation or something
#------------------------------------------------------------------------------#
# Simple function for getting error, accounting for autocorrelation
# def error(r):
#     """
#     Manually calculate the error according to Wilks definition; compare to bse.
#     Proves that the bse parameter is equal to Equation 7.18b in Wilks.
#     """
#     # Use fittedvalues and resid attributes
#     # xfact = (x**2).sum()/(r.nobs*(x-x.mean())**2).sum() # for constant, a
#     # x = r.model.exog[:,-1] # exogeneous variable
#     # xfact = 1/((x-x.mean())**2).sum()
#     se2 = (r.resid.values**2).sum()/(r.nobs-2) # Eq 7.9, Wilks
#     xfact2 = 12/(r.nobs**3-r.nobs) # Eq 3, Thompson et. al 2015
#     sigma = (se2**0.5)*(xfact2**0.5)
#     print('Provided: %.10f, manually calculated: %.10f' % (r.bse.x1, sigma))
#     return

#------------------------------------------------------------------------------#
# Function for performing regression, along with significance
# on regression coefficient or something
# Below uses statsmodels, but slower for small vectors
# X = sm.add_constant(T.index.values) # use the new index values
# r = sm.OLS(T[name,region_get], X).fit() # regress column [name] against X
# scale = ((r.nobs-2)/(r.nobs*((1-autocorr)/(1+autocorr))-2))**0.5
# L[name,region] = [unit*r.params.x1, scale*unit*r.bse.x1]
# Below uses np.corrcoef, but is slower than pearsonr
# autocorr = np.corrcoef(r.resid.values[1:], r.resid.values[:-1])[0,1] # outputs 2by2 matrix
# def regress(x, y, unit=10, sigma=False, ignorenan=False):
#     """
#     Gets regression results with fastest methods possible.
#     See the IPython noteobok for comparisons.
#     """
#     # NaN check
#     x, y = x.squeeze(), y.squeeze()
#     if y.ndim>1:
#         raise ValueError('y is not 1-dimensional.')
#     nan = np.isnan(y)
#     if nan[0] and ignorenan:
#         # If we are getting continuous trends, don't bother computing this one, because
#         # we will also get trends from x[1:], x[2:], etc.
#         return np.nan, np.nan, np.nan, np.nan
#     if nan.any():
#         # Filter out NaNs, get regression on what remains
#         x, y = x[~nan], y[~nan]
#     if x.size<5:
#         # Cannot get decent sigma estimate, in this case
#         return np.nan, np.nan, np.nan, np.nan
#     # Regress, and get sigma if requested
#     # First value is estimated change K through record, second is K/unit
#     if sigma:
#         p, V = np.polyfit(x, y, deg=1, cov=True)
#         resid = y - (x*p[0] + p[1]) # very fast step; don't worry about this one
#         autocorr, _ = st.pearsonr(resid[1:], resid[:-1])
#         scale = (x.size-2)/(x.size*((1-autocorr)/(1+autocorr))-2)
#             # scale factor from Thompson et. al, 2015, Quantifying role of internal variability...
#         stderr = np.sqrt(V[0,0]*scale)
#         return (x[-1]-x[0])*p[0], unit*p[0], (x[-1]-x[0])*stderr, unit*stderr
#     else:
#         p = np.polyfit(x, y, deg=1)
#         return (x[-1]-x[0])*p[0], unit*p[0]

#-------------------------------------------------------------------------------
# TODO: Finish harvesting this original function written for the Objective
# Analysis assignment. Consider deleting it.
#-------------------------------------------------------------------------------
# def spectral(fnm, nm, data, norm=True, win=501,
#             freq_scale=1, scale='days',
#             xlog=True, ylog=False, mc='k',
#             xticks=None,
#             rcolors=('C3','C6'), pcolors=('C0','C1'), alpha=0.99, marker=None,
#             xlim=None, ylim=None, # optional override
#             linewidth=1.5,
#             red_discrete=True, red_contin=True, manual=True, welch=True):
#     '''
#     Spectral transform function. Needs work; was copied from OA
#     assignment and right now just plots a bunch of stuff.
#     '''
#     # Iniital stuff
#     N = len(data)
#     pm = int((win-1)/2)
#     fig, a = plt.subplots(figsize=(15,5))
#
#     # Confidence intervals
#     dof_num = 1.2*2*2*(N/(win/2))
#     trans = 0.5
#     exp = 1
#     F99 = stats.f.ppf(1-(1-alpha)**exp,dof_num,1000)
#     F01 = stats.f.ppf((1-alpha)**exp,dof_num,1000)
#     print('F stats:',F01,F99)
#     rho = np.corrcoef(data[1:],data[:-1])[0,1]
#     kr = np.arange(0,win//2+1)
#     fr = freq_scale*kr/win
#
#     def xtrim(f):
#         if xlim is None:
#             return np.ones(f.size, dtype=bool)
#         else:
#             return ((f>=xlim[0]) & (f<=xlim[-1]))
#
#     # Power spectra
#     if manual:
#         label = 'power spectrum'
#         if welch: label = 'manual method'
#         # Now, manual method with proper overlapping etc.
#         if False:
#             data = data[:(N//pm)*pm]
#         loc = np.linspace(pm,N-pm-1,2*int(np.round(N/win))).round().astype(int) # sample loctaions
#         han = np.hanning(win)
#         han = han/han.sum()
#         phi = np.empty((len(loc),win//2))
#         for i,l in enumerate(loc):
#             pm = int((win-1)/2)
#             C = np.fft.fft(han*signal.detrend(data[l-pm:l+pm+1]))
#             phii = np.abs(C)**2/2
#             phii = 2*phii[1:win//2+1]
#             phi[i,:] = phii
#         phi = phi.mean(axis=0)
#         print('phi sum:',phi.sum())
#         f = np.fft.fftfreq(win)[1:win//2+1]*freq_scale
#         if norm: phi = phi/phi.sum()
#         f, phi = f[xtrim(f)], phi[xtrim(f)] # trim
#         a.plot(f, phi, label=label,
#                mec=mc, mfc=mc, mew=linewidth,
#                marker=marker, color=pcolors[0], linewidth=linewidth)
#         if xlim is None: xlim = ((f*freq_scale).min(), (f*freq_scale).max())
#         if ylim is None: ylim = ((phi.min()*0.95, phi.max()*1.05))
#
#     if welch:
#         label = 'power spectrum'
#         if manual: label = 'welch method'
#         # Welch
#         fw, phi_w = signal.welch(data, nperseg=win, detrend='linear', window='hanning', scaling='spectrum',
#                               return_onesided=False)
#         fw, phi_w = fw[1:win//2+1]*freq_scale, phi_w[1:win//2+1]
#         if norm: phi_w = phi_w/phi_w.sum()
#         fw, phi_w = fw[xtrim(fw)], phi_w[xtrim(fw)] # trim
#         print('phiw sum:',phi_w.sum())
#         a.plot(fw, phi_w, label=label,
#               mec=mc, mfc=mc, mew=linewidth,
#                marker=marker, color=pcolors[-1], linewidth=linewidth)
#         if xlim is None: xlim = ((fw).min(), (fw).max())
#         if ylim is None: ylim = (phi_w.min()*0.95, phi_w.max()*1.05)
#
#     # Best fit red noise spectrum
#     if red_discrete:
#         print('Autocorrelation',rho)
#         phi_r1 = (1-rho**2)/(1+rho**2-2*rho*np.cos(kr*np.pi/(win//2)))
#         print('phi_r1 sum:',phi_r1.sum())
#         if norm: phi_r1 = phi_r1/phi_r1.sum()
#         frp, phi_r1 = fr[xtrim(fr)], phi_r1[xtrim(fr)]
#         a.plot(fr[xtrim(fr)], phi_r1, label=r'red noise, $\rho(\Delta t)$',
#                marker=None, color=rcolors[0], linewidth=linewidth)
#         a.plot(frp, phi_r1*F99, linestyle='--',
#                marker=None, alpha=trans, color=rcolors[0], linewidth=linewidth)
#
#     # Alternate best fit
#     if red_contin:
#         Te = -1/np.log(rho)
#         omega = (kr/win)*np.pi*2
#         phi_r2 = 2*Te/(1+(Te**2)*(omega**2))
#         print('phi_r2 sum:',phi_r2.sum())
#         if norm: phi_r2 = phi_r2/phi_r2.sum()
#         frp, phi_r2 = fr[xtrim(fr)], phi_r2[xtrim(fr)]
#         a.plot(frp, phi_r2, label=r'red noise, $T_e$',
#                marker=None, color=rcolors[1], linewidth=linewidth)
#         a.plot(frp, phi_r2*F99, linestyle='--',
#                marker=None, alpha=trans, color=rcolors[-1], linewidth=linewidth)
#     # Variance
#     print('true variance:',data.std()**2)
#     # Figure formatting
#     a.legend()
#     if ylog:
#         a.set_yscale('log')
#     if xlog:
#         a.set_xscale('log')
#     a.set_title('%s power spectrum' % nm)
#     a.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%5.3g'))
#     a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
#     if xticks is None: xticks = a.get_xticks()
#     my.format(a, xlabel=('frequency (%s${}^{-1}$)' % scale), ylabel='proportion variance explained',
#              xlim=xlim, ylim=ylim, xticks=xticks)
#     suffix = 'pdf'
#     fig.savefig('a5_' + fnm + '.' + suffix, format=suffix, dpi='figure')
#     plt.show()
