def datetime2year(dt):
    """
    Calculate the year fraction from datetime object.
    Input...
        dt (positional):
    """
    # ...pretty straightforward
    thisyear_start = datetime(year=dt.year, month=1, day=1)
    nextyear_start = datetime(year=dt.year+1, month=1, day=1)
    year_part = dt - thisyear_start
    year_length = nextyear_start - thisyear_start
    yearnum  = dt.year + year_part/year_length # we can divide timedelta objects
    return yearnum
