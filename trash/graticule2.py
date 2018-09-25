def graticule(*vectors):
    """
    Return graticule from each vector; used, for example, to make pcolor/pc.
    Assumes evenly spaced data; don't test for this, because of float issues.
    """
    outvectors = tuple()
    for v in vectors:
        delta = v[1]-v[0]
        outvectors += (np.arange(v[0]-delta/2, v[-1]+delta, delta),)
        # outvectors += (np.concatenate((v[:1]-delta/2, (v[1:]+v[:-1])/2, v[-1:]+delta/2)),)
    return outvectors # don't want to just add delta from 0; can get issues

