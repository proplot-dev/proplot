def graticule(lon, lat=None, globe=True):
    """
    Gets graticule for full globe of longitude, latitude data.
    Input...
        lat (positional, required): latitude ndarray vector
        lon (positional, required): longitude ndarray vector
        ...or single arg, and do this for arbitrary vector; works in general
    Output...
        latb: ndarray vector of latitude graticule/mesh/grid boundaries
        lonb: ndarray vector of latitude graticule/mesh/grid boundaries
    DECIDED THIS WAS SILLY; REALLY JUST ACCOMADATES GAUSSIAN GRIDS (WHICH I NEVER USE, 
    AND SHOULD INTERPOLATE TO LATLON ANYWAY), and the WEIRD ERA-INTERIM CONVENTION
    WHERE GRID CELL CENTERS ARE AT +/-90 (WHICH I CAN PROBABLY JUST IGNORE IN
    THE REGION MEANS, BECAUSE THEY ARE ALREADY EXTRA-SMALL CELLS AND CONTRIBUTE
    ON ORDER OF .01 (5deg cells; 2.5deg polar cell; just sum up cosines) TO GLOBE AVERAGE,
    SO DOESN'T MATTER IF IGNORED DUE TO CENTERS AT +/-90... SHOULD PROBABLY JUST DOWNLOAD GAUSSIAN AND
    RE-SAMPLE MYSELF... or ignore it; don't think it should effect basemaps too much... and that
    contribution above is literally worst case scenario, more likely have 1-2 deg resolution and tiny
    contribution, won't even show up in regional averages... maybe just deal with it manually when
    you have era-interim data... i think best case it FIGURE OUT HOW TO DOWNLOAD RAW DATA, AND 
    THEN DOWNSAMPLE IT MYSELF, BECAUSE THEIR METHOD IS WEIRD. GET RAW GAUSSIAN DATA.
    """
    # Enforce monotonic
    if not (lon==np.sort(lon)).all():
        raise ValueError('Longitudes are not monotonic.')
    # Fix lons (or generic vector)
    lodiff = lon[1]-lon[0]
    hidiff = lon[-1]-lon[-2]
    lon = np.concatenate(
            (np.array([lon[0]-lodiff/2]), (lon[1:]+lon[:-1])/2, np.array([lon[-1]+hidiff/2]))
            )

    # Enforce monotonic
    if not (lat==np.sort(lat)).all():
        raise ValueError('Latitudes are not monotonic.')
    # If grid 'center' reported at pole, publishers really mean from the lower 
    # graticule 'up to' the pole; so modify the reported grid center
    if lat[-1]==90:
        lat[-1] = 90-(90-lat[-2])/4
        latbf = 90
    else:
        latbf = lat[-1] + (lat[-1]-lat[-2])/2
    if lat[0]==-90:
        lat[0] = -90.+(lat[1]+90.)/4
        latbi = -90
    else:
        latbi = lat[0] - (lat[1]-lat[0])/2
    lat = np.append(np.append(latbi, (lat[1:]+lat[:-1])/2), latbf)
    # Warnings
    if (lon[0] % 360) != (lon[-1] % 360) and globe:
        print('Warning: output longitude mesh is not circular.')

    # Return stuff
    if lat is None:
        return lon
    else:
        return lon, lat
