def complete_years(dt, startmonth=1):
    """
    THIS THING SUCKS, DON'T USE IT; ALSO VECTORIZE IS A CONVENIENCE FUNCTION, NOT 
    AN EFFICIENCY FUNCTION; MEANT TO BE USED A BILLION TIMES; IMPLEMENTATION IS
    ACTUALLY JUST A PYTHON FOR LOOP; READ DOCUMENTATION.
    Input...
        dt (positional, required): an iterable of datetime objects
        startmonth (kw, optional): the starting month characterizing a "complete"
            year (1 is default, but 12 is also common)
    Returns...
     1) filtered datetimes in original structure, or if structure not indexable
        with a vector, filtered datetime in an ndarray
     2) boolean ndarray vector that can be used to filter out complete years
        given some ndarray of ***monthly*** observation datetimes.
    Test lines...
        import numpy as np
        import pandas as pd
        import myfuncs as my
        from datetime import datetime
        l = [datetime(1970,i,1) for i in range(1,13)]
        d = np.array(l)
        df = pd.DataFrame(np.arange(12),index=d)
        my.completeyears(l)
        my.completeyears(d)
        my.completeyears(df.index)
    """
    # if np.datetime64, or a derived class (e.g. pandas Index with Datetime),
    # should be able to just pull .month from array itself to build array of months
    endmonth = ((startmonth - 2) % 12) + 1
    try:
        # try direct retrieval
        m_idstart = np.where(dt.month==startmonth)[0]
        m_idend = np.where(dt.month==endmonth)[0]
        print('Object is np.datetime64 or derivative.')
    except:
        # try indirect retrieval
        try:
            v = np.vectorize(lambda x: x.month)
            m_idstart = np.where(v(dt)==startmonth)[0]
            m_idend = np.where(v(dt)==endmonth)[0]
            print('Object is ndarray of datetime.datetime objects.')
        except:
            raise TypeError('Invalid input type: %s, and/or subtype: %s'
                    % (type(dt), type(dt[0])))
        pass
    # check if have even one full year
    if m_idstart.size==0 or m_idend.size==0:
        raise ValueError('Input datetime vector does not span at least one full year '
                + 'defined by startmonth %d.' % startmonth)
    # filter out; return from first instance of startmonth to last instance of lastmonth
    filt = np.arange(m_idstart[0], m_idend[-1]+1)
    try:
        # Try simple indexing
        return dt[filt], filt
    except:
        # Cast into array first (since filt is ndarray)
        return np.array(dt)[filt], filt
