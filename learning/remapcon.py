def remapcon(infile, lon, lat, outfile=None):
    """
    Creates formatted file suitable for input as "gridfile" into a CDO interpolation
    scheme, from the vectors of longitude/latitude centers.
    """
    # Check input
    loninc, latinc = lon[1]-lon[0], lat[1]-lat[0]
    
    # Run CDO method
    gridpath = os.path.join(os.path.dirname(infile), 'tmp.grid') # dirname empty, if no slashes
    outfile = os.path.splitext(infile)[0] + '_downsampled.nc' if outfile is None else outfile
        # if dirname is empty, does nothing
    with open(gridpath, 'w+') as grid:
            # the w+ makes new file, if does not exist
            # also w+ overwrites existing content; a+ only appends
        grid.write('\n'.join((
            'gridtype = lonlat',
            'xsize    = %d' % lon.size,
            'ysize    = %d' % lat.size,
            'xfirst   = %.5f' % lon[0], # gives .1 arcsecond precision
            'yfirst   = %.5f' % lat[0],
            'xinc     = %.5f' % loninc,
            'yinc     = %.5f' % latinc,
            ))) # combines strings with '\n' newlines
    subprocess.call((
        'cdo', 
        '-C', # colorized
        '--no_warnings',
        'remapycon,' + gridpath, 
        infile, # input file
        outfile, # output file
        ))
        # use remapycon instead of remapcon; have had issues with latter
    os.remove(gridpath) # always remove
        
    # Return
    print('Name of resampled file: %s.' % outfile)
    return outfile
