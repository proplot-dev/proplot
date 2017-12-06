# ...pre-process options; remove spaces? bad; e.g. days since 1970-1-1 00:00:00
# find "base" with yyyy-m(m)-d(d) or yyyy/m(m)/d(d) search options; each
# month and day can have just one digit
baseregex = re.search(r'\d{4}[-/]\d{1,2}[-/]\d{1,2}', units)
if baseregex is not None: 
    base = datetime.strptime(
            baseregex[0].replace('/','').replace('-',''), '%Y%m%d')
else:
    if verbose: print(errmsg,'time unit string (%s) has ambiguous base date.' % units)
    break
time_nums = M.time
# ...now adjust by converting datenum to timedelta
if 'minutes' in units:
    # make vectorized func for elementwise op -- sometimes need to do this 
    # in numpy, for example when another module doesn't play nice with numpy
    time_nums = time_nums.astype(np.int64)*60 # otherwise get SILENT overflow -- scary!
    time_nums = np.array([timedelta(seconds=int(n)) for n in time_nums.flat]) #np.vectorize(lambda x: timedelta(seconds=int(x)))(time_nums)
elif 'hours' in units:
    # same as above
    time_nums = time_nums.astype(np.int64)*3600 # otherwise get SILENT overflow -- scary!
    time_nums = np.array([timedelta(seconds=int(n)) for n in time_nums.flat])
elif 'days' in units:
    # same as above, but arg is days, not seconds
    time_nums = np.array([timedelta(days=int(n)) for n in time_nums.flat])
else:
    if verbose: print(errmsg,'time unit string (%s) has ambiguous increment unit.' % units)
    break
# ...and add_offset timedelta objects to base time
try: base + time_nums
except OverflowError: # sometimes, datasets with null time or monthly climatology
        # will have dateyear '0', which is below datetime.min == (1,1,1,0,0,0)
    base += timedelta(days=366) # since 'year 0' is leap year
M.time = base + time_nums # will throw error again if something else is wrong
# ...ignore leap years, if requested (models often have no leap years)
if leapyearfix:
    leapyears = np.array([y for y in range(M.time.min().year, M.time.max().year+1) 
        if (y % 4)==0 and ((y % 100)!=0 or (y % 400)==0)]) # numpy min/max work for any datatype allowing comparisons
    leaptimes = np.array([datetime(y, 2, 28, 23, 59, 59) for y in leapyears.flat])
    for l in leaptimes.flat:
        M.time[M.time>l] += timedelta(days=1) # successively add days to these
# ...standardize days for monthly, daily, etc.
if monthly: M.time = np.array([datetime(t.year, t.month, 1, 0) for t in M.time])
if daily: M.time = np.array([datetime(t.year, t.month, t.day, 0) for t in M.time])
if hourly: M.time = np.array([datetime(t.year, t.month, t.day, t.hour) for t in M.time])

