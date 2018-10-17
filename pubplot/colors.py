#!/usr/bin/env python3
"""
File containing tools for converting between colorspaces and various
other color manipulations.
Already provided by matplotlib:
  to_rgb (converts hex-strings, name, rgba to an rgb tuple)
  rgb_to_hsv
  hsv_to_rgb
New utilities:
  hcl_to_rgb
  rgb_to_hcl
  hsluv_to_rgb
  rgb_to_hsluv
"""
import numpy as np
import matplotlib.colors as mcolors

def hcl_to_rgb(c, gamma=3, normalize=False): # gamma is default
    """
    Convert color to HCL space.
    See wiki page: https://en.wikipedia.org/wiki/HCL_color_space
    """
    rgb = mcolors.to_rgb(c) # convert colorish object (e.g. name, hex-code, rgba) to rgb
    alpha = (min(rgb)/max(rgb))/100 # intermediary
    q = np.exp(alpha*gamma) # intermediary
    r, g, b = rgb # expand out, easier
    h = np.arctan2(g-b, r-g)
    h = 2*np.pi+h if h<0 else h # make positive
    h = h*180/np.pi if not normalize else h/(2*np.pi) # normalize to 0-1 possibly
    c = q*(abs(r-g) + abs(g-b) + abs(b-r))/3
    l = (q*max(rgb) + (1-q)*min(rgb))/2
    return (h,c,l)

def shade(color, value=1, saturation=1):
    """
    Modify a color.
    """
    if isinstance(color,str): # will recognize names and hex strings
        color = mcolors.to_rgb(color)
    if any(v>1 for v in color):
        color = (np.array(color)/255).tolist()
    color = mcolors.rgb_to_hsv(color)
    color[2] = max(min(color[2]*value,1),0) # lighten/darken?
    color[1] = max(min(color[1]*saturation,1),0) # saturate/pastelify?
    color = mcolors.hsv_to_rgb(color)
    return mcolors.to_hex(color) # hex codes are nice, mmkay

