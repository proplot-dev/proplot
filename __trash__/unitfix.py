# Probably dangerous, as these are **really really** not standardized, strings 
# can be totally different...better to know a priori the variable types
# Setup
units_old = M.units
scale_factor, add_offset = 1, 0
# ...interpret (run through dict) and apply
flag = False
if M.units.lower() in ('c','degc','deg_c','degreec','degree_c','celsius'):
    pass
elif M.units.lower() in ('k','degk','deg_k','degreek','degree_k','kelvin'):
    M.units, add_offset = 'C', -273.15
elif M.units=='dam':
    pass
elif M.units == 'm':
    M.units, scale_factor = 'dam', 0.1
elif M.units.lower() in ('mb','hpa'):
    pass
elif M.units.lower() == 'pa':
    M.units, scale_factor = 'hPa', 0.01
else:
    flag = True
if flag:
    if verbose: print('Warning: Unknown unit string "%s"' % M.units)
else:
    DATA = DATA*scale_factor + add_offset
    if verbose: print('Data scaled from units "%s" to units "%s" with scale_factor %d and add_offset %d.' 
            % (units_old, M.units, scale_factor, add_offset))
