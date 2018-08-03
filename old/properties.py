def properties(lon, lat):
    """
    Returns properties of graticule/grid underling lon/lat.
    Should also be valid for common Gaussian spectral model-results grids.
    """
    # First, guess cell widths and edges
    dlon1, dlon2, dlat1, dlat2 = lon[1]-lon[0], lon[-1]-lon[-2], lat[1]-lat[0], lat[-1]-lat[-2]
    latb = np.concatenate((lat[:1]-dlat1/2, (lat[1:]+lat[:-1])/2, lat[-1:]+dlat2/2))
    lonb = np.concatenate((lon[:1]-dlon1/2, (lon[1:]+lon[:-1])/2, lon[-1:]+dlon2/2))
    # latb = np.concatenate((lat-dlat/2, lat[-1:]+dlat/2))
    # lonb = np.concatenate((lon-dlon/2, lon[-1:]+dlon/2))
    
    # Cell centers
    latc, lonc = lat.copy(), lon.copy()

    # Corrections
    # ...switch
    has90 = True if lat[-1]==90 else False
    hasm90 = True if lat[0]==-90 else False # need these switches later
    # ...use corrections for dumb grids with 'centers' at poles
    if hasm90:
        latc[0], latb[0] = -90+dlat1/4, -90
    if has90:
        latc[-1], latb[-1] = 90-dlat2/4, 90
    # ...corrected grid widths (cells half as tall near pole)
    # dlat = np.ones(lat.size)*dlat
    # dlon = np.concatenate((np.array([dlon1]), lon[1:]-lon[:-1], np.array([dlon2])))
    # dlat = np.concatenate((np.array([dlat1]), lat[1:]-lat[:-1], np.array([dlat2])))
    dlon = lonb[1:]-lonb[:-1]
    dlat = latb[1:]-latb[:-1]
    if hasm90:
        dlat[0] /= 2
    if has90:
        dlat[-1] /= 2
    
    # Theta/phi coordinates
    phic, phib, dphi = latc*np.pi/180, latb*np.pi/180, dlat*np.pi/180
    thetac, thetab, dtheta = lonc*np.pi/180, lonb*np.pi/180, dlon*np.pi/180

    # Area weights (function of latitude only)
    # ...includes the latitude correction
    areas = dphi[None,:]*dtheta[:,None]*np.cos(phic[None,:])*(c.a**2)
    # areas = dphi*dtheta*np.cos(phic)*(c.a**2)
    #     # close approximation to area; cosine is extremely accurate
    # areas = areas[None,:]
    #     # make lon by lat, so broadcasting rules can apply
    
    # Return dictionary of properties
    # return Dict(
    return Properties(
            latc=latc, latb=latb, dlat=dlat,
            phic=phic, phib=phib, dphi=dphi,
            lonc=lonc, lonb=lonb, dlon=dlon,
            thetac=thetac, thetab=thetab, dtheta=dtheta,
            areas=areas,
            )

