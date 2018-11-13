#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Imports
# Note colormaps are *callable*, will just deliver the corresponding color, easy.
# NOTE:
# * Colormaps generated from HCL space (and cmOcean ones) are indeed perfectly
#   perceptually uniform, but this still looks bad sometimes -- usually we
#   *want* to focus on the *extremes*, so want to weight colors more heavily
#   on the brighters/whiter part of the map! That's what the ColdHot map does,
#   it's what most of the ColorBrewer maps do, and it's what ColorWizard does.
# * By default extremes are stored at end of *lookup table*, not as
#   separate RGBA values (so look under cmap._lut, indexes cmap._i_over and
#   cmap._i_under). You can verify that your cmap is using most extreme values
#   by comparing high-resolution one to low-resolution one.
#------------------------------------------------------------------------------#
# Here's some useful info on colorspaces
# https://en.wikipedia.org/wiki/HSL_and_HSV
# http://www.hclwizard.org/color-scheme/
# http://www.hsluv.org/comparison/ compares lch, hsluv (scaled lch), and hpluv (truncated lch)
# Info on the CIE conventions
# https://en.wikipedia.org/wiki/CIE_1931_color_space
# https://en.wikipedia.org/wiki/CIELUV
# https://en.wikipedia.org/wiki/CIELAB_color_space
# And some useful tools for creating colormaps and cycles
# https://nrlmry.navy.mil/TC.html
# http://help.mail.colostate.edu/tt_o365_imap.aspx
# http://schumacher.atmos.colostate.edu/resources/archivewx.php
# https://coolors.co/
# http://tristen.ca/hcl-picker/#/hlc/12/0.99/C6F67D/0B2026
# http://gka.github.io/palettes/#diverging|c0=darkred,deeppink,lightyellow|c1=lightyellow,lightgreen,teal|steps=13|bez0=1|bez1=1|coL0=1|coL1=1
# https://flowingdata.com/tag/color/
# http://tools.medialab.sciences-po.fr/iwanthue/index.php
# https://learntocodewith.me/posts/color-palette-tools/
#------------------------------------------------------------------------------
import time
import os
import re
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
try:
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a) # noqa
from functools import wraps
from matplotlib import rcParams
from cycler import cycler
from glob import glob
from . import colormath
from . import utils
from .utils import _fill
_data = f'{os.path.dirname(__file__)}' # or parent, but that makes pip install distribution hard
# Define some new palettes
# Note the default listed colormaps
cmap_cycles = ['Set1', 'Set2', 'Set3', 'Set4', 'Set5']
list_cycles = {
    # default matplotlib v2
    'default':      ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    # copied from stylesheets; stylesheets just add color themese from every possible tool, not already present as a colormap
    '538':          ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'],
    'ggplot':       ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
    # the default seaborn ones, they are variations on each other
    'colorblind':   ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9'],
    'colorblind10': ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"], # versions with more colors
    # 'deep':         ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'], # similar to colorblind
    # 'deep10':       ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"],
    # 'muted':        ['#4878CF', '#6ACC65', '#D65F5F', '#B47CC7', '#C4AD66', '#77BEDB'], # similar to colorblind
    # 'muted10':      ["#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4", "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2"],
    # 'bright':       ["#023EFF", "#1AC938", "#E8000B", "#8B2BE2", "#FFC400", "#00D7FF"], # similar to colorblind
    # 'bright10':     ["#023EFF", "#FF7C00", "#1AC938", "#E8000B", "#8B2BE2", "#9F4800", "#F14CC1", "#A3A3A3", "#FFC400", "#00D7FF"],
    # from the website
    'flatui':       ["#3498db", "#e74c3c", "#95a5a6", "#34495e", "#2ecc71", "#9b59b6"],
    # created with online tools
    # add to this!
    # see: http://tools.medialab.sciences-po.fr/iwanthue/index.php
    'cinematic':    [(51,92,103), (158,42,43), (255,243,176), (224,159,62), (84,11,14)],
    'cool':         ["#6C464F", "#9E768F", "#9FA4C4", "#B3CDD1", "#C7F0BD"],
    'sugar':        ["#007EA7", "#B4654A", "#80CED7", "#B3CDD1", "#003249"],
    'vibrant':      ["#007EA7", "#D81159", "#B3CDD1", "#FFBC42", "#0496FF"],
    'office':       ["#252323", "#70798C", "#DAD2BC", "#F5F1ED", "#A99985"],
    'industrial':   ["#38302E", "#6F6866", "#788585", "#BABF95", "#CCDAD1"],
    'tropical':     ["#0D3B66", "#F95738", "#F4D35E", "#FAF0CA", "#EE964B"],
    'intersection': ["#2B4162", "#FA9F42", "#E0E0E2", "#A21817", "#0B6E4F"],
    'field':        ["#23395B", "#D81E5B", "#FFFD98", "#B9E3C6", "#59C9A5"],
    # finally, add the Open Colors ones
    # actually this looks dumb
    # **{'cycle'+str(i): [color+str(i) for color in ('blue', 'red', 'yellow',
    #     'cyan', 'pink', 'teal', 'indigo', 'orange', 'grape', 'lime', 'violet',
    #     'green')] for i in range(10)},
    }
# Aliases
_scale = [359, 99, 99] # for some reason 100 luminance equals 0 luminance!!!
_space_aliases = {
    'rgb':   'rgb',
    'hsv':   'hsv',
    'hpl':   'hpl',
    'hpluv': 'hpl',
    'hsl':   'hsl',
    'hsluv': 'hsl',
    'hcl':   'hcl',
    'lch':   'hcl',
    }
# Names of builtin colormaps
_categories_default = { # initialize as empty lists
    # We keep these ones
    'Matplotlib Originals':
        ['viridis', 'plasma', 'inferno', 'magma', 'twilight', 'twilight_shifted'],
    'SkyPlot Sequential':
        ['Bog', 'Forest', 'Sea', 'Pale', 'Sunrise', 'Sunset', 'Vibrant'], # empty at first, fill automatically
    'SkyPlot Diverging':
        ['ColdHot', 'DryWet', 'Water'],
    'cmOcean Sequential':
        ['Gray', 'Oxy', 'Thermal', 'Haline', 'Ice', 'Dense',
        'Deep', 'Algae', 'Tempo', 'Speed', 'Matter', 'Turbid',
        'Amp', 'Solar', 'Phase', 'Phase_shifted'],
    'cmOcean Diverging':
        ['Balance', 'Curl', 'Delta'],
    'ColorBrewer2.0 Sequential':
        ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
    'ColorBrewer2.0 Diverging':
        ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral'],
    'Other':
        ['cubehelix', 'bwr'],
        # ['cubehelix', 'rainbow', 'bwr'],
    # These ones will be deleted
    'Alt Sequential':
        sorted(['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'coolwarm', 'seismic', # diverging ones
        'afmhot', 'gist_heat', 'copper']),
    'Alt Rainbow':
        sorted(['multi', 'cubehelix', 'cividis']),
    'Alt Diverging':
        sorted(['coolwarm', 'bwr', 'seismic']),
    'Miscellaneous':
        sorted(['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'brg', 'hsv', 'hot', 'rainbow',
        'gist_rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])}
# Categories to ignore/*delete* from dictionary because they suck donkey balls
_categories_ignore = ['Alt Diverging', 'Alt Sequential', 'Alt Rainbow', 'Miscellaneous']

#------------------------------------------------------------------------------#
# More generalized utility for retrieving colors
#------------------------------------------------------------------------------#
def get_space(space):
    """
    Verify requested colorspace is valid.
    """
    space = _space_aliases.get(space, None)
    if space is None:
        raise ValueError(f'Unknown colorspace "{space}".')
    return space

def to_rgb(color, space='rgb'):
    """
    Generalization of mcolors.to_rgb to translate color tuple
    from any colorspace to rgb. Also will convert color strings to tuple.
    """
    # First the RGB input
    # NOTE: Need isinstance here because strings stored in numpy arrays
    # are actually subclasses thereof!
    if isinstance(color, str):
        try:
            color = mcolors.to_rgb(color) # ensure is valid color
        except Exception:
            raise ValueError(f'Invalid RGBA argument {color}. Registered colors are: {", ".join(mcolors._colors_full_map.keys())}.')
    elif space=='rgb':
        color = color[:3] # trim alpha
        if any(c>1 for c in color):
            color = [c/255 for c in color] # scale to within 0-1
    # Next the perceptually uniform versions
    elif space=='hsv':
        color = colormath.hsl_to_rgb(*color)
    elif space=='hpl':
        color = colormath.hpluv_to_rgb(*color)
    elif space=='hsl':
        color = colormath.hsluv_to_rgb(*color)
    elif space=='hcl':
        color = colormath.hcl_to_rgb(*color)
    elif space=='rgb':
        color = color[:3] # trim alpha
        if any(c>1 for c in color):
            color = [c/255 for c in color] # scale to within 0-1
    else:
        raise ValueError('Invalid RGB value.')
    return color

def to_xyz(color, space):
    """
    Inverse of above, translate to some colorspace.
    """
    # Run tuple conversions
    color = mcolors.to_rgb(color) # convert string and/or trim alpha channel
    if space=='hsv':
        color = colormath.rgb_to_hsl(*color) # rgb_to_hsv would also work
    elif space=='hpl':
        color = colormath.rgb_to_hpluv(*color)
    elif space=='hsl':
        color = colormath.rgb_to_hsluv(*color)
    elif space=='hcl':
        color = colormath.rgb_to_hcl(*color)
    elif space=='rgb':
        color = color # do nothing
    else:
        raise ValueError(f'Invalid colorspace {space}.')
    return color

def add_alpha(color):
    """
    Ensures presence of alpha channel.
    """
    if not utils.isvector(color):
        raise ValueError('Input must be color tuple.')
    if len(color)==3:
        color = [*color, 1.0]
    elif len(color)==4:
        color = [*color] # copy, and put into list
    else:
        raise ValueError(f'Tuple length must be 3 or 4, got {len(color)}.')
    return color

def get_channel_value(color, channel, space='hsl'):
    """
    Gets hue, saturation, or luminance channel value from registered
    string color name.

    Arguments
    ---------
        color :
            scalar numeric ranging from 0-1, or string color name, optionally
            with offset specified as '+x' or '-x' at the end of the string for
            arbitrary float x.
        channel :
            channel number or name (e.g., 0, 1, 2, 'h', 's', 'l')
    """
    # Interpret channel
    channel_idxs = {'hue': 0, 'saturation': 1, 'chroma': 1, 'luminance': 2,
                    'alpha': 3, 'h': 0, 's': 1, 'c': 1, 'l': 2}
    channel = channel_idxs.get(channel, channel)
    if callable(color) or not isinstance(color, str):
        return color
    if channel not in (0,1,2):
        raise ValueError('Channel must be in [0,1,2].')
    # Interpret string
    offset = 0
    regex = '([-+]\S*)$' # user can optionally offset from color; don't filter to just numbers, want to raise our own error if user messes up
    match = re.search(regex, color)
    if match:
        try:
            offset = float(match.group(0))
        except ValueError:
            raise ValueError(f'Invalid channel identifier "{color}".')
        color = color[:match.start()]
    return offset + to_xyz(to_rgb(color), space)[channel]

#------------------------------------------------------------------------------#
# Generalized colormap/cycle constructors
#------------------------------------------------------------------------------#
def colormap(*args, extend='both',
        left=None, right=None, x=None, # optionally truncate color range by these indices
        ratios=1, reverse=False, gamma=None, gamma1=None, gamma2=None,
        name=None, register=False, save=False, N=None, **kwargs):
    """
    Convenience function for generating colormaps in a variety of ways.
    The 'extend' property will be used to resample LinearSegmentedColormap
    if we don't intend to use both out-of-bounds colors; otherwise we lose
    the strongest colors at either end of the colormap.

    You can still use extend='neither' in colormap() call with extend='both'
    in contour or colorbar call, just means that colors at ends of the main
    region will be same as out-of-bounds colors.

    Notes on Resampling
    -------------------
    From answer: see https://stackoverflow.com/q/48613920/4970632
    This resampling method is awful! All it does is reduce the
    lookup table size -- what ends up happening under the hood is matplotlib
    tries to *evenly* draw N-1 ('min'/'max') or N-2 ('neither') colors from
    a lookup table with N colors, which means it simply *skips over* 1 or
    2 colors in the middle of the lookup table, which will cause visual jumps!

    Segment data is completely divorced from the number of levels; can
    have many high-res segments with colormap N very small.
    """
    # Turns out pcolormesh makes QuadMesh, which itself is a Collection,
    # which itself gets colors when calling draw() using update_scalarmappable(),
    # which itself uses to_rgba() to get facecolors, which itself is an inherited
    # ScalarMappable method that simply calls the colormap with numbers. Anyway
    # the issue *has* to be with pcolor, because when giving pcolor an actual
    # instance, no longer does that thing where final levels equal extensions.
    # Since collection API does nothing to underlying data or cmap, must be
    # something done by pcolormesh function.
    cmaps = []
    name = name or 'custom' # must have name, mcolors utilities expect this
    N_hires = 256
    if len(args)==0:
        raise ValueError('Function requires at least 1 positional arg.')
    for cmap in args:
        # Retrieve Colormap instance
        if isinstance(cmap,str) and cmap in mcm.cmap_d:
            cmap = mcm.cmap_d[cmap]
        if isinstance(cmap, mcolors.Colormap):
            # Allow gamma override, otherwise do nothing
            if isinstance(cmap, PerceptuallyUniformColormap):
                if gamma1 or gamma2:
                    segmentdata = cmap._segmentdata.copy()
                    if gamma1:
                        segmentdata['gamma1'] = gamma1
                    if gamma2:
                        segmentdata['gamma2'] = gamma2
                    cmap = type(cmap)(cmap.name, segmentdata, space=cmap.space, mask=cmap.mask)
            elif isinstance(cmap, mcolors.LinearSegmentedColormap):
                if gamma:
                    cmap._gamma = gamma
                    cmap._init()
        elif isinstance(cmap, dict):
            # Dictionary of hue/sat/luminance values or 2-tuples representing linear transition
            cmap = PerceptuallyUniformColormap.from_hsl(name, N=N_hires, **cmap)
        elif not isinstance(cmap, str):
            # List of colors
            cmap = mcolors.ListedColormap(cmap, name=name, **kwargs)
        else:
            # Monochrome colormap based from input color (i.e. single hue)
            light = True # by default
            regex = '([0-9].)$'
            match = re.search(regex, cmap) # declare options with _[flags]
            cmap = re.sub(regex, '', cmap) # remove options
            fade = kwargs.pop('fade',90) if not match else match.group(1) # default fade to 90 luminance
            # Build colormap
            cmap = to_rgb(cmap) # to ensure is hex code/registered color
            cmap = monochrome_cmap(cmap, fade, name=name, N=N_hires, **kwargs)
        cmaps += [cmap]
    # Now merge the result of this arbitrary user input
    # Since we are merging cmaps, potentially *many* color transitions; use big number by default
    if len(cmaps)>1:
        N_merge = N_hires*len(cmaps)
        cmap = merge_cmaps(*cmaps, name=name, ratios=ratios, N=N_merge)

    # Reverse
    if reverse:
        cmap = cmap.reversed()

    # Optionally clip edges or resample map.
    try:
        left, right = x
    except TypeError:
        pass
    if isinstance(cmap, mcolors.ListedColormap):
        slicer = None
        if left is not None or right is not None:
            slicer = slice(left,right)
        elif N is not None:
            slicer = slice(None,N)
        # Just sample indices for listed maps
        if slicer:
            slicer = slice(left,right)
            try:
                cmap = mcolors.ListedColormap(cmap.colors[slicer])
            except Exception:
                raise ValueError(f'Invalid indices {slicer} for listed colormap.')
    elif left is not None or right is not None:
        # Trickier for segment data maps
        # First get segmentdata and parse input
        olddata = cmap._segmentdata
        newdata = {}
        if left is None:
            left = 0
        if right is None:
            right = 1
        # Next resample the segmentdata arrays
        for key,xyy in olddata.items():
            if key in ('gamma1', 'gamma2', 'space'):
                newdata[key] = xyy
                continue
            xyy = np.array(xyy)
            x = xyy[:,0]
            xleft, = np.where(x>left)
            xright, = np.where(x<right)
            if len(xleft)==0:
                raise ValueError(f'Invalid x minimum {left}.')
            if len(xright)==0:
                raise ValueError(f'Invalid x maximum {right}.')
            l, r = xleft[0], xright[-1]
            newxyy = xyy[l:r+1,:].copy()
            if l>0:
                xl = xyy[l-1,1:] + (left - x[l-1])*(xyy[l,1:] - xyy[l-1,1:])/(x[l] - x[l-1])
                newxyy = np.concatenate(([[left, *xl]], newxyy), axis=0)
            if r<len(x)-1:
                xr = xyy[r,1:] + (right - x[r])*(xyy[r+1,1:] - xyy[r,1:])/(x[r+1] - x[r])
                newxyy = np.concatenate((newxyy, [[right, *xr]]), axis=0)
            newxyy[:,0] = (newxyy[:,0] - left)/(right - left)
            newdata[key] = newxyy
        # And finally rebuild map
        cmap = type(cmap)(cmap.name, newdata)
    if isinstance(cmap, mcolors.LinearSegmentedColormap) and N is not None:
        # Perform a crude resampling of the data, i.e. just generate a
        # low-resolution lookup table instead
        # NOTE: All this does is create a new colormap with *attribute* N levels,
        # for which '_lut' attribute has not been generated yet.
        offset = {'neither':-1, 'max':0, 'min':0, 'both':1}
        if extend not in offset:
            raise ValueError(f'Unknown extend option {extend}.')
        cmap = cmap._resample(N - offset[extend]) # see mcm.get_cmap source

    # Optionally register a colormap
    if name and register:
        print(name, 'Registering')
        if name.lower() in [cat_cmap.lower() for cat,cat_cmaps in _categories_default.items()
                    for cat_cmap in cat_cmaps if 'SkyPlot' not in cat]:
            print(f'Warning: Overwriting existing colormap "{name}".')
            # raise ValueError(f'Builtin colormap "{name}" already exists. Choose a different name.')
        elif name in mcm.cmap_d:
            pass # no warning necessary
            # print(f'Warning: Overwriting existing colormap "{name}".')
        mcm.cmap_d[name] = cmap
        mcm.cmap_d[name+'_r'] = cmap.reversed()
        if re.search('[A-Z]',name):
            mcm.cmap_d[name.lower()] = cmap
            mcm.cmap_d[name.lower()+'_r'] = cmap.reversed()
        # print(f'Registered name {name}.') # not necessary

    # Optionally save colormap to disk
    if name and save:
        print('Save cmap!')
        # Save segment data directly
        basename = f'{cmap.name}.npy'
        filename = f'{_data}/cmaps/{basename}'
        np.save(filename, dict(cmap._segmentdata, space=cmap.space))
        print(f'Saved colormap to "{basename}".')
        # Save list of hex colors
        # basename = f'{cmap.name}.hex'
        # with open(filename, 'w') as h: # overwrites if exists; otherwise us 'a'
        #     h.write(','.join(mcolors.to_hex(cmap(i)) for i in np.linspace(0,1,cmap.N)))
    return cmap

def colors(*args, vmin=0, vmax=1):
    """
    Convenience function to draw colors from arbitrary ListedColormap or
    LinearSegmentedColormap. Use vmin/vmax to scale your samples.
    """
    samples = 10
    if utils.isnumber(args[-1]):
        args, samples = args[:-1], args[-1]
    elif len(args)>1:
        args = [args] # presumably send a list of colors
    if len(args)==0:
        raise ValueError('Function requires at least 1 positional arg.')
    cmap = colormap(*args) # the cmap object itself
    if isinstance(cmap, mcolors.ListedColormap):
        # Just get the colors
        colors = cmap.colors
    elif isinstance(cmap, mcolors.LinearSegmentedColormap): # or subclass
        # Employ ***more flexible*** version of get_cmap() method, which does this:
        # LinearSegmentedColormap(self.name, self._segmentdata, lutsize)
        if utils.isnumber(samples):
            # samples = np.linspace(0, 1-1/nsample, nsample) # from 'centers'
            samples = np.linspace(0, 1, samples) # from edge to edge
        elif utils.isvector(samples):
            samples = np.array(samples)
        else:
            raise ValueError(f'Invalid samples "{samples}". If you\'re building '
                    'a colormap on-the-fly, input must be [*args, '
                    'samples] where *args are passed to the colormap() constructor '
                    'and "samples" is either the number of samples desired '
                    'or a vector of colormap samples within [0,1].')
        colors = cmap((samples-vmin)/(vmax-vmin))
    else:
        raise ValueError(f'Colormap returned weird object type: {type(cmap)}.')
    return colors

# def colors()
#     """
#     Simple alias.
#     """

class PerceptuallyUniformColormap(mcolors.LinearSegmentedColormap):
    """
    Generate LinearSegmentedColormap in perceptually uniform colorspace --
    i.e. either HSLuv, HCL, or HPLuv. Adds handy feature where *channel
    value for string-name color is looked up*.

    Example
    -------
    dict(hue        = [[0, 'red', 'red'], [1, 'blue', 'blue']],
         saturation = [[0, 1, 1], [1, 1, 1]],
         luminance  = [[0, 1, 1], [1, 0.2, 0.2]])
    """
    def __init__(self, name, segmentdata,
            space='hsl', gamma1=None, gamma2=None,
            mask=False, **kwargs):
        """
        Initialize with dictionary of values. Note that hues should lie in
        range [0,360], saturation/luminance in range [0,100].

        Arguments
        ---------
            mask :
                Whether to mask out-of-range colors as black, or just clip
                the RGB values (distinct from colormap clipping the extremes).
            gamma1 :
                Raise the line used to transition from a low chroma value (x=0)
                to a higher chroma value (x=1) by this power, like HCLWizard.
            gamma2 :
                Raise the line used to transition from a high luminance value (x=0)
                to a lower luminance value (x=1) by this power, like HCLWizard.
        Why change the direction of transition depending on which value is
        bigger? Because makes it much easier to e.g. weight the center of
        a diverging colormap.
        """
        # Attributes
        # NOTE: Don't allow power scaling for hue because that would be weird.
        # Idea is want to allow skewing so dark/saturated colors are
        # more isolated/have greater intensity.
        # NOTE: We add gammas to the segmentdata dictionary so it can be
        # pickled into .npy file
        space = get_space(space)
        if 'gamma' in kwargs:
            raise ValueError('Standard gamma scaling disabled. Use gamma1 or gamma2 instead.')
        segmentdata['gamma1'] = _fill(gamma1, _fill(segmentdata.get('gamma1', None), 1.0))
        segmentdata['gamma2'] = _fill(gamma2, _fill(segmentdata.get('gamma2', None), 1.0))
        self.space = space
        self.mask  = mask
        # First sanitize the segmentdata by converting color strings to their
        # corresponding channel values
        keys   = {*segmentdata.keys()}
        target = {'hue', 'saturation', 'luminance', 'gamma1', 'gamma2'}
        if keys != target and keys != {*target, 'alpha'}:
            raise ValueError('Invalid segmentdata dictionary.')
        for key,array in segmentdata.items():
            # Allow specification of channels using registered string color names
            if 'gamma' in key:
                continue
            if callable(array):
                continue
            for i,xyy in enumerate(array):
                xyy = list(xyy) # make copy!
                for j,y in enumerate(xyy[1:]): # modify the y values
                    xyy[j+1] = get_channel_value(y, key, space)
                segmentdata[key][i] = xyy
        # Initialize
        # NOTE: Our gamma1 and gamma2 scaling is just fancy per-channel
        # gamma scaling, so disable the standard version.
        super().__init__(name, segmentdata, gamma=1.0, **kwargs)

    def reversed(self, name=None):
        """
        Reverse colormap.
        """
        if name is None:
            name = self.name + '_r'
        def factory(dat):
            def func_r(x):
                return dat(1.0 - x)
            return func_r
        data_r = {}
        for key,xyy in self._segmentdata.items():
            if key in ('gamma1', 'gamma2', 'space'):
                if 'gamma' in key: # optional per-segment gamma
                    xyy = np.atleast_1d(xyy)[::-1]
                data_r[key] = xyy
                continue
            elif callable(xyy):
                data_r[key] = factory(xyy)
            else:
                data_r[key] = [[1.0 - x, y1, y0] for x, y0, y1 in reversed(xyy)]
        return PerceptuallyUniformColormap(name, data_r, space=self.space)

    def _init(self):
        """
        As with LinearSegmentedColormap, but convert each value
        in the lookup table from 'input' to RGB.
        """
        # First generate the lookup table
        channels = ('hue','saturation','luminance')
        reverse = (False, False, True) # gamma weights *low chroma* and *high luminance*
        gammas = (1.0, self._segmentdata['gamma1'], self._segmentdata['gamma2'])
        self._lut = np.ones((self.N+3, 4), float) # fill
        for i,(channel,gamma,reverse) in enumerate(zip(channels, gammas, reverse)):
            self._lut[:-3,i] = make_mapping_array(self.N, self._segmentdata[channel], channel, gamma, reverse)
        if 'alpha' in self._segmentdata:
            self._lut[:-3,3] = make_mapping_array(self.N, self._segmentdata['alpha'], 'alpha')
        self._lut[:-3,0] %= 359
        # Make hues circular, and set extremes
        self._isinit = True
        self._set_extremes() # generally just used end values in segmentdata
        # Now convert values to RGBA, and clip colors
        for i in range(self.N+3):
            self._lut[i,:3] = to_rgb(self._lut[i,:3], self.space)
        self._lut[:,:3] = clip_colors(self._lut[:,:3], self.mask)

    def _resample(self, N):
        """
        Return a new color map with *N* entries.
        """
        self.N = N # that easy
        self._i_under = self.N
        self._i_over = self.N + 1
        self._i_bad = self.N + 2
        self._init()
        return self

    @staticmethod
    def from_hsl(name,
            h=0, s=99, l=[99, 20], c=None, a=None,
            hue=None, saturation=None, luminance=None, chroma=None, alpha=None,
            ratios=None, reverse=False, **kwargs):
        """
        Make linear segmented colormap by specifying channel values.
        """
        # Build dictionary, easy peasy
        h = _fill(hue, h)
        s = _fill(chroma, _fill(c, _fill(saturation, s)))
        l = _fill(luminance, l)
        a = _fill(alpha, _fill(a, 1.0))
        cs = ['hue', 'saturation', 'luminance', 'alpha']
        channels = [h, s, l, a]
        cdict = {}
        for c,channel in zip(cs,channels):
            cdict[c] = make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

    @staticmethod
    def from_list(name, color_list,
            ratios=None, reverse=False,
            **kwargs):
        """
        Make linear segmented colormap from list of color tuples. The values
        in a tuple can be strings, in which case that corresponding color-name
        channel value is deduced.

        Optional
        --------
            ratios : simple way to specify x-coordinates for listed color
                transitions -- bigger number is slower transition, smaller
                number is faster transition.
            space : colorspace of hue-saturation-luminance style input
                color tuples.
        """
        # Dictionary
        cdict = {}
        channels = [*zip(*color_list)]
        if len(channels) not in (3,4):
            raise ValueError(f'Bad color list: {color_list}')
        cs = ['hue', 'saturation', 'luminance']
        if len(channels)==4:
            cs += ['alpha']
        else:
            cdict['alpha'] = 1.0 # dummy function that always returns 1.0
        # Build data arrays
        for c,channel in zip(cs,channels):
            cdict[c] = make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

def make_segmentdata_array(values, ratios=None, reverse=False, **kwargs):
    """
    Construct a list of linear segments for an individual channel.
    This was made so that user can input e.g. a callable function for
    one channel, but request linear interpolation for another one.
    """
    # Handle function handles
    if callable(values):
        if reverse:
            values = lambda x: values(1-x)
        return values # just return the callable
    values = np.atleast_1d(values)
    if len(values)==1:
        value = values[0]
        return [(0, value, value), (1, value, value)] # just return a constant transition

    # Get x coordinates
    if not np.iterable(values):
        raise TypeError('Colors must be iterable.')
    if ratios is not None:
        xvals = np.atleast_1d(ratios) # could be ratios=1, i.e. dummy
        if len(xvals) != len(values) - 1:
            raise ValueError(f'Got {len(values)} values, but {len(ratios)} ratios.')
        xvals = np.concatenate(([0], np.cumsum(xvals)))
        xvals = xvals/np.max(xvals) # normalize to 0-1
    else:
        xvals = np.linspace(0,1,len(values))

    # Build vector
    array = []
    slicer = slice(None,None,-1) if reverse else slice(None)
    for x,value in zip(xvals,values[slicer]):
        array.append((x, value, value))
    return array

def make_mapping_array(N, data, channel, gamma=1.0, reverse=False):
    """
    Mostly a copy of matplotlib version, with a few modifications:
        * Permit ranges outside of 0-1 (i.e. allowing HSL values).
        * Allow circular hue gradations.
        * Allow weighting each transition by going from:
            c = c1 + x*(c2 - c1)
          for x in range [0-1], to
            c = c1 + (x**gamma)*(c2 - c1)
    """
    # Optionally allow for ***callable*** instead of linearly interpolating
    # between line segments
    gammas = np.atleast_1d(gamma)
    if (gammas < 0.01).any() or (gammas > 10).any():
        raise ValueError('Gamma can only be in range [0.01,10].')
    if callable(data):
        if len(gammas)>1:
            raise ValueError('Only one gamma allowed for functional segmentdata.')
        x = np.linspace(0, 1, N)**gamma
        lut = np.array(data(x), dtype=float)
        return lut

    # Get array
    try:
        data = np.array(data)
    except Exception:
        raise TypeError('Data must be convertible to an array.')
    shape = data.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError('Data must be nx3 format.')
    if len(gammas)!=1 and len(gammas)!=shape[0]-1:
        raise ValueError(f'Need {shape[0]-1} gammas for {shape[0]}-level mapping array, but got {len(gamma)}.')
    if len(gammas)==1:
        gammas = np.repeat(gammas, shape[:1])

    # Get indices
    x  = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    if x[0] != 0.0 or x[-1] != 1.0:
        raise ValueError('Data mapping points must start with x=0 and end with x=1.')
    if (np.diff(x) < 0).any():
        raise ValueError('Data mapping points must have x in increasing order.')
    x = x*(N - 1)

    # Get distances from the segmentdata entry to the *left* for each requested
    # level, excluding ends at (0,1), which must exactly match segmentdata ends
    xq = (N - 1)*np.linspace(0, 1, N)
    ind = np.searchsorted(x, xq)[1:-1] # where xq[i] must be inserted so it is larger than x[ind[i]-1] but smaller than x[ind[i]]
    distance = (xq[1:-1] - x[ind - 1])/(x[ind] - x[ind - 1])
    # Scale distances in each segment by input gamma
    _, uind, cind = np.unique(ind, return_index=True, return_counts=True)
    for i,(ui,ci) in enumerate(zip(uind,cind)): # i will range from 0 to N-2
        # Test if 1
        gamma = gammas[i]
        if gamma==1:
            continue
        # By default, weight toward a *lower* channel value (i.e. bigger
        # exponent implies more colors at lower value)
        ir = ((y0[i + 1] - y1[i]) < 0) # by default want to weight toward a *lower* channel value
        if reverse:
            ir = (not ir)
        # print('Left value', y1[i], 'Right value', y0[i + 1], 'Reverse', ir)
        if ir:
            distance[ui:ui + ci] = 1 - (1 - distance[ui:ui + ci])**gamma
        else:
            distance[ui:ui + ci] **= gamma

    # Perform successive linear interpolations all rolled up into one equation
    lut = np.zeros((N,), float)
    lut[1:-1] = distance*(y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0]  = y1[0]
    lut[-1] = y0[-1]
    if channel in ('saturation', 'luminance'):
        lut[:] = np.clip(lut[:], 0, 99)
    return lut

#------------------------------------------------------------------------------#
# Colormap constructors
#------------------------------------------------------------------------------#
def merge_cmaps(*cmaps, name='merged', N=512, ratios=1, **kwargs):
    """
    Merge arbitrary colormaps.
    Arguments
    ---------
        cmaps : 
            List of colormap strings or instances for merging.
        name :
            Name of output colormap.
        N :
            Number of lookup table colors desired for output colormap.
    Notes
    -----
    * Old method had us simply calling the colormap with arrays of fractions.
      This was sloppy, because it just samples locations on the lookup table and
      will therefore degrade the original, smooth, functional transitions.
    * Better method is to combine the _segmentdata arrays and simply scale
      the x coordinates in each (x,y1,y2) channel-tuple according to the ratios.
    * In the case of ListedColormaps, we just combine the colors.
    """
    # Initial
    if len(cmaps)<=1:
        raise ValueError('Need two or more input cmaps.')
    ratios = ratios or 1
    if utils.isscalar(ratios):
        ratios = [1]*len(cmaps)

    # Combine the colors
    cmaps = [colormap(cmap, N=None, **kwargs) for cmap in cmaps] # set N=None to disable resamping
    if all(isinstance(cmap,mcolors.ListedColormap) for cmap in cmaps):
        if not np.all(ratios==1):
            raise ValueError(f'Cannot assign different ratios when mering ListedColormaps.')
        colors = [color for cmap.colors in cmaps for color in cmap.colors]
        cmap = mcolors.ListedColormap(colors, name=name, N=len(colors))

    # Accurate methods for cmaps with continuous/functional transitions
    elif all(isinstance(cmap,mcolors.LinearSegmentedColormap) for cmap in cmaps):
        # Combine the actual segmentdata
        kinds = {type(cmap) for cmap in cmaps}
        if len(kinds)>1:
            raise ValueError(f'Got mixed colormap types.')
        kind = kinds.pop() # colormap kind
        keys = {key for cmap in cmaps for key in cmap._segmentdata.keys()}
        ratios = np.array(ratios)/np.sum(ratios) # so if 4 cmaps, will be 1/4
        x0 = np.concatenate([[0], np.cumsum(ratios)])
        xw = x0[1:] - x0[:-1]

        # Combine the segmentdata, and use the y1/y2 slots at merge points
        # so the transition is immediate (can never interpolate between end
        # colors on the two colormaps)
        segmentdata = {}
        gamma1, gamma2 = [], []
        for key in keys:
            # Combine scalar values
            if key in ('gamma1', 'gamma2'):
                if key not in segmentdata:
                    segmentdata[key] = []
                for cmap in cmaps:
                    segmentdata[key] += [cmap._segmentdata[key]]
                continue
            # Combine xyy data
            datas = []
            test = [callable(cmap._segmentdata[key]) for cmap in cmaps]
            if not all(test) and any(test):
                raise ValueError('Mixed callable and non-callable colormap values.')
            if all(test): # expand range from x-to-w to 0-1
                for x,w,cmap in zip(x0[:-1],xw,cmaps):
                    data = lambda x: data((x - x0)/w) # WARNING: untested!
                    datas.append(data)
                def data(x):
                    idx, = np.where(x<x0)
                    if idx.size==0:
                        i = 0
                    elif idx.size==x0.size:
                        i = x0.size-2
                    else:
                        i = idx[-1]
                    return datas[i](x)
            else:
                for x,w,cmap in zip(x0[:-1],xw,cmaps):
                    data = np.array(cmap._segmentdata[key])
                    data[:,0] = x + w*data[:,0]
                    datas.append(data)
                for i in range(len(datas)-1):
                    datas[i][-1,2] = datas[i+1][0,2]
                    datas[i+1] = datas[i+1][1:,:]
                data = np.concatenate(datas, axis=0)
                data[:,0] = data[:,0]/data[:,0].max(axis=0) # scale to make maximum exactly 1 (avoid floating point errors)
            segmentdata[key] = data

        # Create object
        kwargs = {}
        if kind is PerceptuallyUniformColormap:
            spaces = {cmap.space for cmap in cmaps}
            if len(spaces)>1:
                raise ValueError(f'Trying to merge colormaps with different HSL spaces {repr(spaces)}.')
            kwargs.update({'space':spaces.pop()})
        cmap = kind(name, segmentdata, N=N, **kwargs)
    else:
        raise ValueError('All colormaps should be of the same type (Listed or LinearSegmented).')
    return cmap

def monochrome_cmap(color, fade, reverse=False, space='hsl', name='monochrome', **kwargs):
    """
    Make a sequential colormap that blends from color to near-white.
    Arguments
    ---------
        color :
            Build colormap by varying the luminance of some RGB color while
            keeping its saturation and hue constant.
    Optional
    --------
        reverse : (False)
            Optionally reverse colormap.
        space : ('hsl')
            Colorspace in which we vary luminance.
    """
    # Get colorspace
    space = get_space(space)
    h, s, l = to_xyz(to_rgb(color), space)
    if utils.isnumber(fade): # allow just specifying the luminance channel
        fade = np.clip(fade, 0, 99) # 99 is the max!
        fade = to_rgb((h, 0, fade), space=space)
    _, fs, fl = to_xyz(to_rgb(fade), space)
    fs = s # consider changing this?
    index = slice(None,None,-1) if reverse else slice(None)
    return PerceptuallyUniformColormap.from_hsl(name, h, [s,fs][index], [l,fl][index], space=space, **kwargs)

def clip_colors(colors, mask=True):
    """
    Arguments
    ---------
        colors :
            List of length-3 RGB color tuples.
        mask : (bool)
            Whether to mask out (set to some dark gray color) or clip (limit
            range of each channel to [0,1]) out-of-range RGB channels.
    Notes
    -----
    Could use np.clip (matplotlib.colors uses this under the hood) but want
    to display messages, and anyway premature efficiency is the root of all
    evil, we're manipulating like 1000 colors max here, it's no big deal.
    """
    under = {}
    over = {}
    message = 'Invalid value for' if mask else 'Clipping'
    colors = [list(rgb) for rgb in colors] # so we can overwrite
    gray = 0.2 # RGB values for gray color
    for rgb in colors:
        for i,(c,n) in enumerate(zip(rgb,'rgb')):
            if c<0:
                if mask:
                    rgb[-1] = 0
                    for j in range(3):
                        rgb[j] = gray
                else:
                    rgb[i] = 0
                if n not in under:
                    under[n] = True
                    # print(f'Warning: {message} channel {n} (<0).')
            if c>1:
                if mask:
                    rgb[-1] = 1 # full alpha
                    for j in range(3):
                        rgb[j] = gray
                else:
                    rgb[i] = 1 # keep
                if n not in over:
                    over[n] = True
                    # print(f'Warning: {message} channel {n} (>1).')
    return colors

#------------------------------------------------------------------------------#
# Cycle helper functions
#------------------------------------------------------------------------------#
def set_cycle(cmap, samples=None, rename=False):
    """
    Set the color cycler.
    Arguments
    ---------
        cmap :
            Name of colormap or colormap instance from which we draw list of colors.
        samples :
            Array of values from 0-1 or number indicating number of evenly spaced
            samples from 0-1 from which to draw colormap colors. Will be ignored
            if the colormap is a ListedColormap (interpolation not possible).
    """
    colors = colors(cmap, samples)
    cyl = cycler('color', colors)
    rcParams['axes.prop_cycle'] = cyl
    rcParams['patch.facecolor'] = colors[0]
    if rename:
        rename_colors(cmap)

def rename_colors(cycle='colorblind'):
    """
    Calling this will change how shorthand codes like "b" or "g"
    are interpreted by matplotlib in subsequent plots.
    Arguments
    ---------
        cycle : {deep, muted, pastel, dark, bright, colorblind}
            Named seaborn palette to use as the source of colors.
    """
    seaborn_cycles = ['colorblind', 'deep', 'muted', 'bright']
    if cycle=='reset':
        colors = [(0.0, 0.0, 1.0), (0.0, .50, 0.0), (1.0, 0.0, 0.0), (.75, .75, 0.0),
                  (.75, .75, 0.0), (0.0, .75, .75), (0.0, 0.0, 0.0)]
    elif cycle in seaborn_cycles:
        colors = cycles[cycle] + [(0.1, 0.1, 0.1)]
    else:
        raise ValueError(f'Cannot set colors with color cycle {cycle}.')
    for code, color in zip('bgrmyck', colors):
        rgb = mcolors.colorConverter.to_rgb(color)
        mcolors.colorConverter.colors[code] = rgb
        mcolors.colorConverter.cache[code]  = rgb

#------------------------------------------------------------------------------
# Normalization classes for mapping data to colors (i.e. colormaps)
# WARNING: Many methods in ColorBarBase tests for class membership, crucially
# including _process_values(), which if it doesn't detect BoundaryNorm will
# end up trying to infer boundaries from inverse() method
#------------------------------------------------------------------------------
def Norm(norm, **kwargs):
    """
    Return arbitrary normalizer.

    Notes
    -----
    Tried following the pcolor example (see Colormap notes) to use BoundaryNorm
    for getting segmented levels from colormap. Note this is sort-of built into
    mcolors.Colormap (by generating lookup table with N colors) -- BoundaryNorm
    is an alternative. It seems that using normalizers to discretize colors
    performs *better* than using low-resolution resampling when the colormap
    is complex, but not exactly sure why.
    """
    if isinstance(norm, mcolors.Normalize):
        pass
    elif norm is None:
        # norm = 'linear' # always the default
        norm = None
    elif norm is None:
        norm = mcolors.Normalize() # default is just linear from 0 to 1
    elif type(norm) is not str: # dictionary lookup
        raise ValueError(f'Unknown norm "{norm}".')
    if isinstance(norm, str):
        if norm not in normalizers:
            raise ValueError(f'Unknown normalizer "{norm}". Options are {", ".join(normalizers.keys())}.')
        norm = normalizers[norm](**kwargs)
    return norm

class LinearSegmentedNorm(mcolors.Normalize):
    """
    Description
    -----------
    As in BoundaryNorm case, but instead we linearly *interpolate* colors
    between the provided boundary levels. Exactly analagous to the method
    in LinearSegmentedColormap: perform linear interpolations between
    successive monotonic, but arbitrarily spaced, points.
    """
    def __init__(self, levels, norm=None, midpoint=None, clip=False, ncolors=None, **kwargs):
        # Very simple
        levels = np.atleast_1d(levels)
        if levels.size<=1 or ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError(f'Levels passed to LinearSegmentedNorm must be monotonically increasing.')
        super().__init__(np.nanmin(levels), np.nanmax(levels), clip) # second level superclass
        self._x = levels # alias for boundaries
        # self._pre_norm = norm # e.g. a log-norm
        # if norm:
        #     levels = norm(levels)

    def __call__(self, xq, clip=None):
        # Follow example of make_mapping_array for efficient, vectorized
        # linear interpolation across multiple segments
        # NOTE: normal test puts values at a[i] if a[i-1] < v <= a[i]; for
        # left-most data, satisfy a[0] <= v <= a[1]
        # NOTE: searchsorted gives where xq[i] must be inserted so it is larger
        # than x[ind[i]-1] but smaller than x[ind[i]]
        x = self._x # from arbitrarily spaced monotonic levels
        y = np.linspace(0, 1, x.size) # to linear range 0-1
        xq = np.atleast_1d(xq)
        ind = np.searchsorted(x, xq)
        ind[ind==0] = 1
        distance = (xq - x[ind - 1])/(x[ind] - x[ind - 1])
        yq = distance*(y[ind] - y[ind - 1]) + y[ind - 1]
        return ma.masked_array(yq, np.isnan(xq))

    def inverse(self, yq):
        # Performs inverse operation of __call__
        x = self._x
        y = np.linspace(0, 1, x.size)
        yq = np.atleast_1d(yq)
        ind = np.searchsorted(y, yq)
        ind[ind==0] = 1
        distance = (yq - y[ind - 1])/(y[ind] - y[ind - 1])
        xq = distance*(x[ind] - x[ind - 1]) + x[ind - 1]
        return ma.masked_array(xq, np.isnan(yq))
        # if self._pre_norm:
        #     x = self._pre_norm.inverse(x) # e.g. first scale logarithmically

class BinNorm(mcolors.BoundaryNorm):
    """
    Simple normalizer that *interpolates* from an RGB array at point
    (level_idx/num_levels) along the array, instead of choosing color
    from (transform(level_value)-transform(vmin))/(transform(vmax)-transform(vmin))
    where transform can be linear, logarithmic, etc.

    Note
    ----
    If you are using a diverging colormap with extend='max/min', the center
    will get messed up. But that is very strange usage anyway... so please
    just don't do that :)

    Todo
    ----
    Allow this to accept transforms too, which will help prevent level edges
    from being skewed toward left or right in case of logarithmic/exponential data.

    Example
    -------
    Your levels edges are weirdly spaced [-1000, 100, 0, 100, 1000] or
    even [0, 10, 12, 20, 22], but center "colors" are always at colormap
    coordinates [.2, .4, .6, .8] no matter the spacing; levels just must be monotonic.
    """
    def __init__(self, levels, norm=None, centers=False, clip=False, extend=None, **kwargs):
        # NOTE: Idea is that we bin data into len(levels) discrete x-coordinates,
        # and optionally make out-of-bounds colors the same or different
        # NOTE: Don't need to call parent __init__, this is own implementation
        # Do need it to subclass BoundaryNorm, so ColorbarBase will detect it
        # See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        # Declare boundaries, vmin, vmax in True coordinates
        extend = extend or 'both'
        levels = np.atleast_1d(levels)
        if levels.size<=1 or ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError('Levels passed to Normalize() must be monotonically increasing.')
        if extend not in ('both','min','max','neither'):
            raise ValueError(f'Unknown extend option "{extend}". Choose from "min", "max", "both", "neither".')
        if centers:
            levels = utils.edges(levels)
        N = len(levels)
        # Determine y-bin *centers* desired
        # NOTE: Length of bin centers should be N + 1
        norm = norm or (lambda x: x) # e.g. a logarithmic transform
        offset = {'both':2, 'min':1, 'max':1, 'neither':0}
        resample = {'both':range(N+1), 'neither':[0, *range(N-1), N-2],
                    'min':[*range(N), N-1], 'max':[0, *range(N)]}
        self._x = norm(levels)
        self._y_bins = np.linspace(0, 1, N + offset[extend] - 1)[resample[extend]]
        self._pre_norm = norm
        # Add builtin properties
        self.boundaries = levels # alias read by other functions
        self.vmin = levels[0]
        self.vmax = levels[-1]
        self.clip = clip
        self.N = N

    def __call__(self, xq, clip=None):
        # Follow example of LinearSegmentedNorm, but perform no interpolation,
        # just use searchsorted to bin the data
        # NOTE: The bins vector includes out-of-bounds negative (searchsorted
        # index 0) and out-of-bounds positive (searchsorted index N+1) values
        x = self._pre_norm(self._x)
        xq = np.atleast_1d(xq)
        yq = self._y_bins[np.searchsorted(x, xq)] # which x-bin does each point in xq belong to?
        return ma.masked_array(yq, np.isnan(xq))

    def inverse(self, yq):
        # Not possible
        raise ValueError('BinNorm is not invertible.')

class StretchNorm(mcolors.Normalize):
    """
    Normalizers that 'stretches' and 'compresses' either side of a colormap
    about some midpoint, proceeding exponentially (exp>0) or logarithmically
    (exp<0) down the linear colormap from the center point.

    Notes
    -----
    * Default midpoint is vmin, i.e. we just stretch to the right. For diverging
      colormaps, use midpoint 0.5.
    * Need to update this. Should features be incorporated with
      LinearSegmentedNorm? Should user just use exponential gradation
      functions in the segmentdata instead of using a special normalizer?
    """
    def __init__(self, exp=0, midpoint=None, vmin=None, vmax=None, clip=None):
        # Bigger numbers are too one-sided
        if abs(exp) > 10:
            raise ValueError('Warping scale must be between -10 and 10.')
        super().__init__(vmin, vmax, clip)
        self._midpoint = midpoint
        self._exp = exp

    def _warp(x, exp, exp_max=4):
        # Returns indices stretched so neutral/low values are sampled more heavily
        if exp > 0:
            invert = True
        else:
            invert, exp = False, -exp
        exp = exp*(exp_max/10)
        # Apply function; approaches x=1 as a-->Inf, x=x as a-->0
        if invert: x = 1-x
        value =  (x-1+(np.exp(x)-x)**exp)/(np.e-1)**exp
        if invert:
            value = 1-value # flip on y-axis
        return value

    def __call__(self, value, clip=None):
        # Get middle point in 0-1 coords, and value
        midpoint = self._midpoint or self.vmin
        midpoint_scaled = (midpoint - self.vmin)/(self.vmax - self.vmin)
        value_scaled    = (value - self.vmin)/(self.vmax - self.vmin)
        try: iter(value_scaled)
        except TypeError:
            value_scaled = np.arange(value_scaled)
        value_cmap = ma.empty(value_scaled.size)
        # Get values, accounting for midpoints
        for i,v in enumerate(value_scaled):
            v = np.clip(v, 0, 1)
            if v>=midpoint_scaled:
                block_width = 1 - midpoint_scaled
                value_cmap[i] = (midpoint_scaled +
                    block_width*self._warp((v - midpoint_scaled)/block_width, self._exp)
                    )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled -
                    block_width*self._warp((midpoint_scaled - v)/block_width, self._exp)
                    )
        return value_cmap

#------------------------------------------------------------------------------#
# Register new colormaps; must come before registering the color cycles
# * If leave 'name' empty in register_cmap, name will be taken from the
#   Colormap instance. So do that.
# * Note that **calls to cmap instance do not interpolate values**; this is only
#   done by specifying levels in contourf call, specifying lut in get_cmap,
#   and using LinearSegmentedColormap.from_list with some N.
# * The cmap object itself only **picks colors closest to the "correct" one
#   in a "lookup table**; using lut in get_cmap interpolates lookup table.
#   See LinearSegmentedColormap doc: https://matplotlib.org/api/_as_gen/matplotlib.colors.LinearSegmentedColormap.html#matplotlib.colors.LinearSegmentedColormap
# * If you want to always disable interpolation, use ListedColormap. This type
#   of colormap instance will choose nearest-neighbors when using get_cmap, levels, etc.
#------------------------------------------------------------------------------#
def register_colors(nmax=np.inf, threshold=0.10):
    """
    Register new color names. Will only read first n of these
    colors, since XKCD library is massive (they should be sorted by popularity
    so later ones are no loss).

    Notes
    -----
    * The 'threshold' arg denotes how separated each channel of the HCL converted
        colors must be.
    * This seems like it would be slow, but takes on average 0.03 seconds on
        my macbook, so it's fine.
    """
    # First ***reset*** the colors dictionary
    # Why? We want to add XKCD colors *sorted by popularity* from file, along
    # with crayons dictionary; having registered colors named 'xkcd:color' is
    # annoying and not useful
    translate =  {'b': 'blue', 'g': 'green', 'r': 'red', 'c': 'cyan',
                  'm': 'magenta', 'y': 'yellow', 'k': 'black', 'w': 'white'}
    base1 = mcolors.BASE_COLORS # one-character names
    base2 = {translate[key]:value for key,value in base1.items()} # full names
    mcolors._colors_full_map.clear() # clean out!
    mcolors._colors_full_map.cache.clear() # clean out!
    mcolors._colors_full_map.update(base1)
    mcolors._colors_full_map.update(base2)

    # First register colors and get their HSL values
    space = 'hcl'
    hcls = np.empty((0,3))
    names = []
    categories_list = []
    categories_set = set()
    for file in glob(f'{_data}/colors/*.txt'):
        # Read data
        category, _ = os.path.splitext(os.path.basename(file))
        data = np.genfromtxt(file, delimiter='\t', dtype=str, comments='%', usecols=(0,1)).tolist()
        ncolors = min(len(data),nmax-1)
        hcl = np.empty((ncolors,3))
        # Add categories
        categories_set.add(category)
        custom_colors[category] = {}
        custom_colors_filtered[category] = {}
        # Translate and put in dictionary
        for i,(name,color) in enumerate(data): # is list of name, color tuples
            if i>=nmax: # e.g. for xkcd colors
                break
            hcl[i,:] = to_xyz(color, space=space)
            name = re.sub('grey', 'gray', name)
            name = re.sub('/', ' ', name)
            name = re.sub(r'\bpinky\b', 'pink', name)
            name = re.sub(r'\bgreeny\b', 'green', name)
            name = re.sub("'s", '', name)
            names.append((category, name))
            custom_colors[category][name] = color
        # Concatenate HCL arrays
        hcls = np.concatenate((hcls, hcl), axis=0)

    # Remove colors that are 'too similar' by rounding to the nearest n units
    # WARNING: unique axis argument requires numpy version >=1.13
    # WARNING: evidently it is ***impossible*** to actually delete colors
    # from the custom_colors dictionary (perhaps due to quirk of autoreload,
    # perhaps by some more fundamental python thing), so we instead must create
    # *completely separate* dictionary and add colors from there
    hcls = hcls/np.array(_scale)
    hcls = np.round(hcls/threshold).astype(np.int64)
    _, index, counts = np.unique(hcls, return_index=True, return_counts=True, axis=0) # get unique rows
    deleted = 0
    counts = counts.sum()
    exceptions = ['white', 'black', 'gray', 'red', 'pink', 'grape',
            'violet', 'indigo', 'blue', 'coral', 'tomato red', 'crimson',
            'cyan', 'teal', 'green', 'lime', 'yellow', 'orange']
    exceptions_regex = '^(' + '|'.join(exceptions) + ')[0-9]?$'

    # Add colors to filtered colors
    for i,(category,name) in enumerate(names):
        if not re.match(exceptions_regex, name) and i not in index:
            deleted += 1
        else:
            custom_colors_filtered[category][name] = custom_colors[category][name]
    for category,dictionary in custom_colors_filtered.items():
        mcolors._colors_full_map.update(dictionary)
    # print(f'Started with {len(names)} colors, removed {deleted} insufficiently distinct colors.')

def register_cmaps():
    """
    Register colormaps and cycles in the cmaps directory.
    Note all of those methods simply modify the dictionary mcm.cmap_d.
    """
    # First read from file
    for file in glob(f'{_data}/cmaps/*'):
        # Read table of RGB values
        if not re.search('.(rgb|hex|npy)$', file):
            continue
        name = os.path.basename(file)[:-4]
        # Comment this out to overwrite existing ones
        # if name in mcm.cmap_d: # don't want to re-register every time
        #     continue
        if re.search('.rgb$', file):
            try: cmap = np.loadtxt(file, delimiter=',') # simple
            except:
                print(f'Failed to load {os.path.basename(file)}.')
                continue
            if (cmap>1).any():
                cmap = cmap/255
        # Read list of hex strings
        elif re.search('.hex$', file):
            cmap = [*open(file)] # just a single line
            if len(cmap)==0:
                continue # file is empty
            cmap = cmap[0].strip().split(',') # csv hex strings
            cmap = np.array([mcolors.to_rgb(c) for c in cmap]) # from list of tuples
        # Directly read segmentdata of hex strings
        # Will ensure that HSL colormaps have the 'space' entry
        else:
            segmentdata = np.load(file).item() # unpack 0-D array
            if 'space' in segmentdata:
                space = segmentdata.pop('space')
                cmap = PerceptuallyUniformColormap(name, segmentdata, N=N, space=space)
            else:
                cmap = mcolors.LinearSegmentedColormap(name, segmentdata, N=N)
        # Register as ListedColormap or LinearSegmentedColormap
        if isinstance(cmap, mcolors.Colormap): # i.e. we did not load segmentdata directly
            cmap_r = cmap.reversed()
        else:
            if 'lines' in name.lower():
                cmap   = mcolors.ListedColormap(cmap)
                cmap_r = cmap.reversed() # default name is name+'_r'
                custom_cycles.add(name)
            else:
                N = len(cmap) # simple as that; number of rows of colors
                cmap   = mcolors.LinearSegmentedColormap.from_list(name, cmap, N) # using static method is way easier
                cmap_r = cmap.reversed() # default name is name+'_r'
                custom_cmaps.add(name)
        # Register maps (this is just what register_cmap does)
        mcm.cmap_d[cmap.name]   = cmap
        mcm.cmap_d[cmap_r.name] = cmap_r

    # Fix the builtin rainbow colormaps by switching from Listed to
    # LinearSegmented -- don't know why matplotlib shifts with these as
    # discrete maps by default, dumb.
    for name in _categories_default['Matplotlib Originals']: # initialize as empty lists
        cmap = mcm.cmap_d.get(name,None)
        if cmap and isinstance(cmap, mcolors.ListedColormap):
            mcm.cmap_d[name] = mcolors.LinearSegmentedColormap.from_list(name, cmap.colors)

    # Swap the order of divering colorbrewer maps, direction of color changes
    # is opposite from intuition (red to blue, pink to green, etc.)
    names = []
    if 'RdBu' in mcm.cmap_d: # only do this once! we modified the content of ColorBrewer Diverging
        for name in _categories_default['ColorBrewer2.0 Diverging']:
            # Reverse map and name
            # e.g. RdBu --> BuRd, RdYlBu --> BuYlRd
            # Note default name PuOr is literally backwards...
            cmap   = mcm.cmap_d.get(name, None)
            cmap_r = mcm.cmap_d.get(name + '_r', None)
            if cmap:
                if name not in ('Spectral','PuOr','BrBG'):
                    del mcm.cmap_d[name]
                    del mcm.cmap_d[name + '_r']
                    name = re.sub('^(..)(..)?(..)$', r'\3\2\1', name)
                    mcm.cmap_d[name] = cmap_r
                    mcm.cmap_d[name + '_r'] = cmap
            names += [name]
    _categories_default['ColorBrewer2.0 Diverging'] = names

    # Add shifted versions of cyclic colormaps, and prevent same colors on ends
    for name in ['twilight', 'Phase']:
        cmap = mcm.cmap_d.get(name, None)
        if cmap and isinstance(cmap, mcolors.LinearSegmentedColormap):
            data = cmap._segmentdata
            data_shift = data.copy()
            for key,array in data.items():
                array = np.array(array)
                # Drop an end color
                array = array[1:,:]
                array_shift = array.copy()
                array_shift[:,0] -= 0.5
                array_shift[:,0] %= 1
                array_shift = array_shift[array_shift[:,0].argsort(),:]
                # Normalize x-range
                array[:,0] -= array[:,0].min()
                array[:,0] /= array[:,0].max()
                data[key] = array
                array_shift[:,0] -= array_shift[:,0].min()
                array_shift[:,0] /= array_shift[:,0].max()
                data_shift[key] = array_shift
            # Register shifted version and original
            mcm.cmap_d[name] = mcolors.LinearSegmentedColormap(name, data, cmap.N)
            mcm.cmap_d[name + '_shifted'] = mcolors.LinearSegmentedColormap(name + '_shifted', data_shift, cmap.N)

    # Delete ugly ones
    for category in _categories_ignore:
        for name in _categories_default:
            mcm.cmap_d.pop(name, None)

    # Register names so that they can be invoked ***without capitalization***
    # This always bugged me! Note cannot change dictionary during iteration.
    ignorecase = {}
    ignore = [category for categories in _categories_ignore for category in categories]
    for name,cmap in mcm.cmap_d.items():
        if name in ignore:
            mcm.cmap_d.pop(name, None)
        elif re.search('[A-Z]',name):
            ignorecase[name.lower()] = cmap
    mcm.cmap_d.update(ignorecase)
    for key in ignorecase.keys():
        custom_lower.add(key)

def register_cycles():
    """
    Register cycles defined right here by dictionaries.
    """
    # Simply register them as ListedColormaps
    for name,colors in list_cycles.items():
        mcm.cmap_d[name]        = mcolors.ListedColormap([to_rgb(color) for color in colors])
        mcm.cmap_d[f'{name}_r'] = mcolors.ListedColormap([to_rgb(color) for color in colors[::-1]])

    # Remove some redundant ones
    mcm.cmap_d.pop('tab10', None)
    mcm.cmap_d.pop('tab20', None)
    mcm.cmap_d.pop('Paired', None)
    mcm.cmap_d.pop('Pastel1', None)
    mcm.cmap_d.pop('Pastel2', None)
    mcm.cmap_d.pop('Dark2', None)
    if 'Accent' in mcm.cmap_d:
        mcm.cmap_d.pop('Set1', None)
        mcm.cmap_d['Set1'] = mcm.cmap_d.pop('Accent')
    if 'tab20b' in mcm.cmap_d:
        mcm.cmap_d['Set4'] = mcm.cmap_d.pop('tab20b')
    if 'tab20c' in mcm.cmap_d:
        mcm.cmap_d['Set5'] = mcm.cmap_d.pop('tab20c')

# Register stuff when this module is imported
# The 'cycles' are simply listed colormaps, and the 'cmaps' are the smoothly
# varying LinearSegmentedColormap instances or subclasses thereof
custom_cmaps = set() # track downloaded colormaps
custom_cycles = set() # same, but track cycles
custom_lower = set() # lower-case keys added to dictionary, that we will ignore
custom_colors = {} # downloaded colors categorized by filename
custom_colors_filtered = {} # limit to 'sufficiently unique' color names
register_cmaps()
register_colors()
register_cycles()
print('Registered colors and colormaps.')

# Finally our dictionary of normalizers
# Includes some custom classes, so has to go at end
normalizers = {
    'none':       mcolors.NoNorm,
    'null':       mcolors.NoNorm,
    'step':       BinNorm,
    'segmented':  LinearSegmentedNorm,
    'boundary':   mcolors.BoundaryNorm,
    'log':        mcolors.LogNorm,
    'linear':     mcolors.Normalize,
    'power':      mcolors.PowerNorm,
    'symlog':     mcolors.SymLogNorm,
    }

#------------------------------------------------------------------------------#
# Visualizations
#------------------------------------------------------------------------------#
def color_show(groups=['open', ['crayons','xkcd']], ncols=4, nbreak=12, minsat=0.2):
    """
    Visualize all possible named colors. Wheee!
    Modified from: https://matplotlib.org/examples/color/named_colors.html
    * Special Note: The 'Tableau Colors' are just the *default matplotlib
      color cycle colors*! So don't bother iterating over them.
    """
    # Get colors explicitly defined in _colors_full_map, or the default
    # components of that map (see soure code; is just a dictionary wrapper
    # on some simple lists)
    figs = []
    for group in groups:
        # Get group colors
        group = group or 'open'
        if isinstance(group, str):
            group = [group]
        color_dict = {}
        for name in group:
            # Read colors from current cycler
            if name=='cycle':
                seen = set() # trickery
                cycle_colors = rcParams['axes.prop_cycle'].by_key()['color']
                cycle_colors = [color for color in cycle_colors if not (color in seen or seen.add(color))] # trickery
                color_dict.update({f'C{i}':v for i,v in enumerate(cycle_colors)})
            # Read custom defined colors
            else:
                color_dict.update(custom_colors_filtered[name]) # add category dictionary

        # Group colors together by discrete range of hue, then sort by value
        # For opencolors this is not necessary
        if 'open' in group:
            # Sorted color columns and plot settings
            space = 0.5
            swatch = 1.5
            names = ['gray', 'red', 'pink', 'grape', 'violet', 'indigo', 'blue', 'cyan', 'teal', 'green', 'lime', 'yellow', 'orange']
            nrows, ncols = 10, len(names) # rows and columns
            plot_names = [[name+str(i) for i in range(nrows)] for name in names]
        # For other palettes this is necessary
        else:
            # Get colors in perceptally uniform space
            # Then will group based on hue thresholds
            space = 1
            swatch = 1
            colors_hsl = {key:
                [channel/scale for channel,scale in zip(to_xyz(value, 'hcl'), _scale)]
                for key,value in color_dict.items()}

            # Keep in separate columns
            breakpoints = np.linspace(0,1,nbreak) # group in blocks of 20 hues
            plot_names = [] # initialize
            sat_test = (lambda x: x<minsat) # test saturation for 'grays'
            for n in range(len(breakpoints)):
                # Get 'grays' column
                if n==0:
                    hue_colors = [(name,hsl) for name,hsl in colors_hsl.items()
                                  if sat_test(hsl[1])]
                # Get column for nth color
                else:
                    b1, b2 = breakpoints[n-1], breakpoints[n]
                    hue_test   = ((lambda x: b1<=x<=b2) if b2 is breakpoints[-1]
                                   else (lambda x: b1<=x<b2))
                    hue_colors = [(name,hsl) for name,hsl in colors_hsl.items() if
                            hue_test(hsl[0]) and not sat_test(hsl[1])] # grays have separate category
                # Get indices to build sorted list, then append sorted list
                sorted_index = np.argsort([pair[1][2] for pair in hue_colors])
                plot_names.append([hue_colors[i][0] for i in sorted_index])
            # Concatenate those columns so get nice rectangle
            # nrows = max(len(huelist) for huelist in plot_names) # number of rows
            ncols = nbreak-1
            names = [i for sublist in plot_names for i in sublist]
            plot_names = [[]]
            nrows = len(names)//ncols+1
            for i,name in enumerate(names):
                if ((i + 1) % nrows)==0:
                    plot_names.append([]) # add new empty list
                plot_names[-1].append(name)

        # Create plot by iterating over columns 
        # Easy peasy
        figsize = (8*space*(ncols/4), 5*(nrows/40)) # 5in tall with 40 colors in column
        fig, ax = plt.subplots(figsize=figsize)
        X, Y = fig.get_dpi()*fig.get_size_inches() # size in *dots*; make these axes units
        h, w = Y/(nrows+1), X/ncols # height and width of row/column in *dots*
        for col,huelist in enumerate(plot_names):
            for row,name in enumerate(huelist): # list of colors in hue category
                y = Y - (row * h) - h
                xi_line = w*(col + 0.05)
                xf_line = w*(col + 0.25*swatch)
                xi_text = w*(col + 0.25*swatch + 0.03*swatch)
                print_name = name.split('xkcd:')[-1] # make sure no xkcd:
                ax.text(xi_text, y, print_name, fontsize=h*0.8, ha='left', va='center')
                ax.hlines(y+h*0.1, xi_line, xf_line, color=color_dict[name], lw=h*0.6)
        ax.set_xlim(0,X)
        ax.set_ylim(0,Y)
        ax.set_axis_off()
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)

        # Save figure
        fig.savefig(f'{_data}/colors/colors_{"-".join(group)}.pdf',
                bbox_inches='tight', format='pdf', transparent=False)
        figs += [fig]
    return figs

def cycle_show():
    """
    Show off the different color cycles.
    Wrote this one myself, so it uses the custom API.
    """
    # Get the list of cycles
    # cycles = plt.get_cyclS
    cycles = {**{name:mcm.cmap_d[name].colors for name in cmap_cycles},
              **{name:mcm.cmap_d[name].colors for name in list_cycles.keys()}}
    nrows = len(cycles)//2+len(cycles)%2
    # Create plot
    fig, axs = plt.subplots(figsize=(6,nrows*1.5), ncols=2, nrows=nrows)
    fig.subplots_adjust(
            top=0.98, bottom=0.01, left=0.02, right=0.98,
            hspace=0.25, wspace=0.02)
    axs = [ax for sub in axs for ax in sub]
    state = np.random.RandomState(528)
    for i,(ax,(key,cycle)) in enumerate(zip(axs,cycles.items())):
        lw = 4
        key = key.lower()
        array = state.rand(20,len(cycle)) - 0.5
        array = array[:,:1] + array.cumsum(axis=0) + np.arange(0,len(cycle))
        for j,color in enumerate(cycle):
            l, = ax.plot(array[:,j], lw=5, ls='-', color=color)
            l.set_zorder(10+len(cycle)-j) # make first lines have big zorder
        title = f'{key}: {len(cycle)} colors'
        ax.set_title(title)
        ax.grid(True)
        for axis in 'xy':
            ax.tick_params(axis=axis,
                    which='both', labelbottom=False, labelleft=False,
                    bottom=False, top=False, left=False, right=False)
    if len(cycles)%2==1:
        axs[-1].set_visible(False)
    # Save
    fig.savefig(f'{_data}/colors/cycles.pdf',
            bbox_inches='tight', format='pdf')
    return fig

# def cmap_show(N=31, ignore=['Miscellaneous','Sequential2','Diverging2']):
def cmap_show(N=31):
    """
    Plot all current colormaps, along with their catgories.
    This example comes from the Cookbook on www.scipy.org. According to the
    history, Andrew Straw did the conversion from an old page, but it is
    unclear who the original author is.
    See: http://matplotlib.org/examples/color/colormaps_reference.html
    """
    # Have colormaps separated into categories:
    # NOTE: viridis, cividis, plasma, inferno, and magma are all
    # listed colormaps for some reason
    exceptions = ['viridis','cividis','plasma','inferno','magma']
    track_lowercase = set()
    cmaps_all = [cmap for cmap in mcm.cmap_d.keys() if
            not cmap.endswith('_r')
            and cmap not in custom_lower
            and 'Vega' not in cmap
            and (isinstance(mcm.cmap_d[cmap],mcolors.LinearSegmentedColormap) or cmap in exceptions)]
    cmaps_listed = [cmap for cmap in mcm.cmap_d.keys() if
            not cmap.endswith('_r')
            and cmap not in custom_lower
            and 'Vega' not in cmap
            and (not isinstance(mcm.cmap_d[cmap],mcolors.LinearSegmentedColormap) and cmap not in exceptions)]

    # Detect unknown/manually created colormaps, and filter out
    # colormaps belonging to certain section
    categories    = {cat:cmaps for cat,cmaps in _categories_default.items() if cat not in _categories_ignore}
    cmaps_ignore  = [cmap for cat,cmaps in _categories_default.items() for cmap in cmaps if cat in _categories_ignore]
    cmaps_missing = [cmap for cat,cmaps in categories.items() for cmap in cmaps if cmap not in cmaps_all]
    cmaps_known   = [cmap for cat,cmaps in categories.items() for cmap in cmaps if cmap in cmaps_all]
    cmaps_custom  = [cmap for cmap in cmaps_all if cmap not in cmaps_known and cmap not in cmaps_ignore]

    # Attempt to auto-detect diverging colormaps, just sample the points on either end
    # Do this by simply summing the RGB channels to get HSV brightness
    # l = lambda i: to_xyz(to_rgb(m(i)), 'hcl')[2] # get luminance
    # if (l(0)<l(0.5) and l(1)<l(0.5)): # or (l(0)>l(0.5) and l(1)>l(0.5)):
    # if cmap.lower() in custom_diverging:
    if cmaps_missing:
        print(f'Missing colormaps: {", ".join(cmaps_missing)}')
    if cmaps_ignore:
        print(f'Ignored colormaps: {", ".join(cmaps_ignore)}')
    if cmaps_custom:
        print(f'New colormaps: {", ".join(cmaps_custom)}')

    # Attempt sorting based on hue
    # for cat in ['SkyPlot Sequential', 'cmOcean Sequential', 'ColorBrewer2.0 Sequential']:
    # for cat in ['SkyPlot Sequential', 'ColorBrewer2.0 Sequential']:
    for cat in []:
        hues = [np.mean([to_xyz(to_rgb(color),'hsl')[0]
            for color in mcm.cmap_d[cmap](np.linspace(0.3,1,20))])
            for cmap in categories[cat]]
        categories[cat] = [categories[cat][idx] for idx,cmap in zip(np.argsort(hues),categories[cat])]

    # Array for producing visualization with imshow
    a = np.linspace(0, 1, 257).reshape(1,-1)
    a = np.vstack((a,a))

    # Figure
    extra = 1 # number of axes-widths to allocate for titles
    nmaps = len(cmaps_known) + len(cmaps_custom) + len(categories)*extra
    fig = plt.figure(figsize=(5,0.3*nmaps))
    fig.subplots_adjust(top=0.98, bottom=0.005, left=0.20, right=0.99)

    # Make plot
    ntitles, nplots = 0, 0 # for deciding which axes to plot in
    for cat in categories:
        # Space for title
        ntitles += extra # two axes-widths
        for i,m in enumerate(categories[cat]):
            # Checks
            if i+ntitles+nplots>nmaps:
                break
            if m not in mcm.cmap_d or m not in cmaps_all: # i.e. the expected builtin colormap is missing
                continue
            # Draw, and make axes invisible
            cmap = mcm.get_cmap(m, N) # interpolate
            ax = plt.subplot(nmaps,1,i+ntitles+nplots)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.imshow(a, aspect='auto', cmap=cmap, origin='lower')
            # Category title and colorbar label
            if i==0:
                t = ax.title
                t.set_text(cat) # category name
                t.set_visible(True)
            yl = ax.yaxis.label
            yl.set_text(m) # map name
            yl.set_visible(True)
            yl.set_rotation(0)
            yl.set_ha('right')
            yl.set_va('center')
        # Space for plots
        nplots += len(categories[cat])

    # Save
    filename = f'{_data}/cmaps/colormaps.pdf'
    print(f"Saving figure to: {filename}.")
    fig.savefig(filename, bbox_inches='tight')
    return fig

