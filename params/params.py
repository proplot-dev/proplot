#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#
import os
import numpy as np
import xarray as xr
import scipy.signal as signal
import scipy.stats as stats
import time as T # have some big algorithms here; want to guage how long they take
from . import const

#------------------------------------------------------------------------------#
# Cross-section and zonal-mean parameters
# Perhaps should always rely in XArray objects
# TODO: Put stuff here from timescales experiment.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------
# TODO: Physical meteorological quantities
# This sections is *very limited*, and should consider deleting it
# Should probably be replaced with *climate.py* stuff
#------------------------------------------------------------------------------
def pt(lev, T, axis=2):
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

