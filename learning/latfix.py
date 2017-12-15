def latfix(lat, data=None):
    """
    Fixed lattiudes.
    Copied from block in ncutils.
    This is stupid simple, but lets me be brainless.
    Completely silly, forget this idea.
    """
    if data is not None:
        if data.shape[0]!=lat.size:
            raise ValueError('Data must have latitude on second dimension.')
    # Flip latitudes
    if lat[1]<lat[0]:
       lat = np.flip(lat)
       if data is not None:
           data = np.flip(data,axis=1)
    if data is None:
        return lat
    else:
        return lat, data
