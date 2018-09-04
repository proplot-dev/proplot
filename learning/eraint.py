def eraint(variables, mode, date_range, 
        level_type='sfc', level_range=None, levels=None,
        hour=(0,6,12,18),
        res=0.75, box=None,
        filename='DATA.nc'):
    """
    Retrieves ERA-Interim DATA using the provided API.
    variables: can be either of...
        list/tuple of variable string names
        individual variable string
    * must know MARS id for requested variables; can add to dictionary in code below using
        https://rda.ucar.edu/datasets/ds627.0/docs/era_interim_grib_table.html
    mode: can be any of...
        'synoptic' (6 hourly DATA)
        'monthly' (monthly mean)
        'accumulation' (things like precipitation; has units per day)
    date_range: can be either of...
        yyyy-mm-dd format string, for single day
        length-2 tuple/list of yyyy-mm-dd strings, data will span these
    level_type: can be any of...
        'p' (pressure levels)
        'sfc' (earth surface)
        'theta', 'pt' (potential temperature)
        'pv' (2pvu surface)
    level_range: can be either of...
        length-2 tuple/list of pressure/pt levels; retrieves all available levels between these
        single number, to pick individual level
    levels: can be either of...
        list/tuple of multiple levels
        single number, to pick individual level 
    hour: can be either of...
        list/tuple of integer hours (should be in [0,6,12,18])
        single number, to pick individual hour
    res: desired output resolution; not sure if this is arbitrary, or if ERA-interim only has
        a select few valid resolution options.
    box: can be either of...
        string name for particular region, e.g. "europe" (see documentation)
        the West/South/East/North boundaries (so lower-left corner, upper-right corner), as a length-4 list/tuple
    filename: name of file output
    """
    # Variable id conversion (see: https://rda.ucar.edu/datasets/ds627.0/docs/era_interim_grib_table.html)
    variables_convert = dict(
            t2m='167.128', # 2m temp
            d2m='168.128', # 2m dew point
            sst='34.128', # sst
            msl='151.128', # sea level pressure
            slp='151.128', # same
            z='129.128', # geopotential
            t='130.128', # temp
            u='131.128', # u wind
            v='132.128', # v wind
            w='135.128', # w wind
            q='133.128', # specific humidity
            r='157.128', # relative humidity
            pt='3.128', # potential temp (available on 2pvu surf)
            theta='3.128', # same
            p='54.128', # pressure (availble on pt, 2pvu surfaces)
            pres='54.128', # same
            pv='60.128' # potential vorticity (available on p, pt surfaces)
            )
    if type(variables) is str:
        variables = (variables,)
    for v in variables:
        if v not in variables_convert:
            raise ValueError('MARS id for variable "%s" is unknown (might need to be added to this script).' % v)
    variables = '/'.join(variables_convert[v] for v in variables)

    # Time range conversion
    if type(date_range) is str:
        date_range = (date_range, date_range) # put into list
    date_range = '/to/'.join(date_range)

    # Hour conversion
    try:
        iter(hour)
    except TypeError: # put into list, if singleton
        hour = (hour,)
    hour = '/'.join(str(h).zfill(2) for h in hour)

    # Data mode
    mode_convert = dict(synoptic='oper',monthly='moda',accumulation='mdfa')
    if mode not in mode_convert:
        raise ValueError('Unknown DATA retrieval mode: "%s"' % mode)
    mode = mode_convert[mode]
        # oper is original, moda is monthly mean of daily means, monthly is 
        # monthly mean of synoptic mean (e.g. 0Z), accumulation is monthly mean
        # of things like rain in units *per day*

    # Level type
    type_convert = dict(p='pl',sfc='sfc',theta='pt',pt='pt',pv='pv')
    # type_full = dict(
    #         p=(1,2,3,5,7,10,20,30,50,70,100,125,150,175,200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000),
    #         sfc=None,
    #         theta=(265,270,285,300,315,330,350,370,395,430,475,530,600,700,850),
    #         pt=(265,270,285,300,315,330,350,370,395,430,475,530,600,700,850),
    #         pv=(2000,)
    #         )
    if level_type not in type_convert:
        raise ValueError('Unknown level type: "%s"' % level_type)
    level_type = type_convert[level_type]

    # Levels - allow specific selection, or range from what's available
    # ...convert levels
    try:
        if levels is not None: iter(levels)
    except TypeError:
        levels = (levels,)
    # ...or convert level_range
    try:
        if level_range is not None: iter(level_range)
    except TypeError:
        levels = (level_range,)
    if levels is not None:
        levels = '/'.join(str(l) for l in levels)
    if level_range is not None: # if user specified both, will use level_range
        levels = '/to/'.join(str(l) for l in level_range)
            # MARS will fill the /to/ with increments of one

    # Area - can be specified as pre-defined region (e.g. string 'europe') OR n/s/w/e boundary
    if box is not None and type(box) is not str:
        box = [box[3], box[0], box[2], box[1]] # MARS takes north/west/east/south
        box = '/'.join(str(a) for a in box)

    # Main server instructions
    # ...not really sure what happens in some situations: list so far...
    # 1) evidently if you provide with variable string-name instead of numeric ID, 
    #       MARS will search for correct one; if there is name ambiguity/conflict will throw error
    # 2) on GUI framework, ECMWF only offers a few resolution options, but program seems
    #       to run when requesting custom resolutions like 5deg/5deg
    retrieve = {
        'class':    'ei', # ecmwf classifiction; choose ERA-Interim
        'type':     'an', # type of field; analysis 'an' or forecast 'fc'
        'format':   'netcdf', # want netcdf instead of GRIB
        'stream':   mode, # product monthly, raw, etc.
        'date':     date_range,
        'time':     hour,
        'levtype':  level_type,
        'param':    variables,
        'dataset':  'interim',
        'step':     '0', # number of hours forecast has been run into future from 'time'
        'grid':     str(res) + '/' + str(res),
        'target':   filename, # save location
        'resol':    'av', # prevents truncation before transformation to geo grid
    }

    # Optional/conditional instructions
    if levels is not None:
        retrieve['levelist'] = levels
    if box is not None:
        retrieve['area'] = box
    if mode[0]!='m':
        retrieve['hour'] = hour
    print('Final MARS request: %s' % retrieve)

    # Retrieve DATA with settings
    server = ECMWFDataServer()
    server.retrieve(retrieve)
    return
