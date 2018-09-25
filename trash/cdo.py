def cdo(process, infiles, args=None, outfile=None):
    """
    OLD VERSION; DOESN'T WORK GREAT FOR PIPES.
    Calls CDO process, with arguments specified by args.
    Just a convenience function so I can call CDO methods from python scripts.
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
    # Create command
    infiles = (infiles,) if type(infiles) is str else infiles
    args = '' if args is None else ','+args if type(args) is str else ','+','.join(args)
    outfile = '%s_%s.nc' % (os.path.splitext(infiles[0])[0], process) if outfile is None else outfile # output file
    sub.call((
        'cdo',
        '-C', # colorized
        process + args,  # process with arguments
        *infiles, # input file
        outfile
        ))
    
    # And return outfile
    return outfile
