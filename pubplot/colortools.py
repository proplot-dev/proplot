#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Imports
# Note colormaps are *callable*, will just deliver the corresponding color, easy.
# TODO: Still confused over some issues:
# * Note that by default extremes are stored at end of *lookup table*, not as
#   separate RGBA values (so look under cmap._lut, indexes cmap._i_over and
#   cmap._i_under). You can verify that your cmap is using most extreme values
#   by comparing high-resolution one to low-resolution one.
# * Seems by default that extremes are always a separate color; so far have
#   been using *hack* where we simply *resample the lookup table* to the number
#   of levels minus extremes, and magically extremes work out to the same color.
# * Unsure whether to (a) limit the _segmentdata in a segmented colormap to
#   the desired number with get_cmap() (which calls _resample), (b) create
#   a ListedColormap with fixed colors, or (d) simply create a discrete
#   normalizer that rounds to the nearest color in high-resolution lookup
#   table? Probably should prefer the latter?
# * Unsure of issue with contourf, StepNorm, and colorbar that results in
#   weird offset ticks. Contour function seems to do weird stuff, still has
#   high-resolution LinearSegmentedColormap, but distinct contour levels. Maybe
#   you shouldn't pass a discrete normalizer since contour takes care of this
#   task itself.
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
#------------------------------------------------------------------------------#
import os
import re
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from functools import wraps
from matplotlib import rcParams
from cycler import cycler
from glob import glob
from . import colormath
from . import utils
# Define some new palettes
cycles = {
    # default matplotlib v2
    'default':      ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    # copied from stylesheets; stylesheets just add color themese from every possible tool, not already present as a colormap
    'ggplot':       ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
    'bmh':          ['#348ABD', '#A60628', '#7A68A6', '#467821', '#D55E00', '#CC79A7', '#56B4E9', '#009E73', '#F0E442', '#0072B2'],
    'solarized':    ['#268BD2', '#2AA198', '#859900', '#B58900', '#CB4B16', '#DC322F', '#D33682', '#6C71C4'],
    '538':          ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'],
    'seaborn':      ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'],
    'pastel':       ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6'],
    'colorblind':   ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9'],
    'deep':         ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'], # similar to colorblind
    'muted':        ['#4878CF', '#6ACC65', '#D65F5F', '#B47CC7', '#C4AD66', '#77BEDB'], # similar to colorblind
    'bright':       ["#023EFF", "#1AC938", "#E8000B", "#8B2BE2", "#FFC400", "#00D7FF"], # similar to colorblind
    'colorblind10': ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"], # versions with more colors
    'deep10':       ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"],
    'muted10':      ["#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4", "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2"],
    'bright10':     ["#023EFF", "#FF7C00", "#1AC938", "#E8000B", "#8B2BE2", "#9F4800", "#F14CC1", "#A3A3A3", "#FFC400", "#00D7FF"],
    # created with online tools
    'cinematic1':   [(51,92,103), (255,243,176), (224,159,62), (158,42,43), (84,11,14)],
    'cinematic2':   [(1,116,152), (231,80,0), (123,65,75), (197,207,255), (241,255,47)],
    }
seaborn_cycles = ['colorblind', 'deep', 'muted', 'bright']
# Note the default listed colormaps
cmap_cycles = ['Pastel1', 'Pastel2', 'Paired',
    'Accent', 'Dark2', 'Set1', 'Set2', 'Set3',
    'tab10', 'tab20', 'tab20b', 'tab20c']
# Finally some names
space_aliases = {
    'rgb':   'rgb',
    'hsv':   'hsv',
    'hpl':   'hpl',
    'hpluv': 'hpl',
    'hsl':   'hsl',
    'hsluv': 'hsl',
    'hcl':   'hcl',
    'lch':   'hcl',
    }
# Scale
space_scales = {
    'rgb': (1,1,1,1),
    'hsv': (1,1,1,1),
    'hsl': (359,99,99,1),
    'hpl': (359,99,99,1),
    'hcl': (359,99,99,1),
    }
# Aliases
channel_idxs = {'hue': 0, 'saturation': 1, 'chroma': 1, 'luminance': 2, 'alpha': 3,
                'h':   0, 's':          1, 'c':      1, 'l':         2}
# Names of builtin colormaps
categories_default = { # initialize as empty lists
    'Builtin Rainbow':
        sorted(['viridis', 'plasma', 'inferno', 'magma']),
    'Rainbow Alt':
        sorted(['multi', 'cubehelix', 'cividis']),
    'Custom Sequential': [], # should come first
    'Custom Diverging': [], # should come first
    # 'ColorBrewer': [], # actually *every single colorbrewer map* is implemented already by default, just different names
    'Builtin Diverging':
        sorted(['piyg', 'prgn', 'brbg', 'puor', 'rdgy', 'rdbu', 'rdylbu', 'rdylgn', 'spectral']),
    'Builtin Sequential':
        sorted(['greys', 'purples', 'blues', 'greens', 'oranges', 'reds',
        'ylorbr', 'ylorrd', 'orrd', 'purd', 'rdpu', 'bupu',
        'gnbu', 'pubu', 'ylgnbu', 'pubugn', 'bugn', 'ylgn']),
    'Sequential Alt':
        sorted(['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'wistia',
        'coolwarm', 'bwr', 'seismic', # diverging ones
        'afmhot', 'gist_heat', 'copper']),
    'Diverging Alt': sorted(['coolwarm', 'bwr', 'seismic']),
    'Miscellaneous':
        sorted(['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'cmrmap', 'brg', 'hsv', 'hot',
        'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])}

#------------------------------------------------------------------------------#
# More generalized utility for retrieving colors
#------------------------------------------------------------------------------#
def get_space(space):
    """
    Verify requested colorspace is valid.
    """
    space = space_aliases.get(space, None)
    if space is None:
        raise ValueError(f'Unknown colorspace "{space}".')
    return space

def get_scale(space):
    """
    Get scales.
    """
    return space_scales[get_space(space)]

def to_rgb(color, space='rgb'):
    """
    Generalization of mcolors.to_rgb to translate color tuple
    from any colorspace to rgb. Also will convert color strings to tuple.
    """
    # First the RGB input
    if type(color) is str:
        color = mcolors.to_rgb(color) # ensure is valid color
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
        color = colormath.rgb_to_hsl(*color)
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

def get_channel(color, channel, space='hsl', scale=True):
    """
    Get hue, saturation, or luminance from string color name. If input
    is not string, will just scale it appropriately for the colorspace.
    Arguments
    ---------
        color : scalar numeric ranging from 0-1, or string color name, optionally
            with offset specified as '+x' or '-x' at the end of the string for
            arbitrary float x.
        channel : channel number (can be 0, 1, or 2).
    Optional
    --------
        scale : whether input tuples need to be scaled to full ranges
            i.e. multiply by 99, 99, 359 for most HSl colorspaces.
    """
    # Interpret channel
    channel = channel_idxs.get(channel, channel)
    if channel not in (0,1,2,3):
        raise ValueError('Channel must be in [0,1,2].')
    # Bail out
    if utils.isvector(color):
        raise TypeError('Input should be string or scalar number.')
    if type(color) is not str:
        scales = get_scale(space) if scale else (1,1,1,1)
        return color*scales[channel]
    if channel==3:
        raise ValueError(f'Cannot specify alpha channel with color string.')
    # Interpret string
    scales = get_scale(space)
    offset = 0
    regex = '([-+]\S*)$' # user can optionally offset from color; don't filter to just numbers, want to raise our own error if user messes up
    match = re.search(regex, color)
    if match:
        try:
            offset = float(match.group(0))
        except ValueError as err:
            raise type(err)(f'Invalid channel identifier "{color}".')
        color = color[:match.start()]
    return scales[channel]*offset + to_xyz(to_rgb(color), space)[channel]

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
def register_colors(nmax=256):
    """
    Register new color names. Will only read first n of these
    colors, since XKCD library is massive (they should be sorted by popularity
    so later ones are no loss).
    """
    for file in glob(f'{os.path.dirname(__file__)}/colors/*.txt'):
        category, _ = os.path.splitext(os.path.basename(file))
        data = np.genfromtxt(file, delimiter='\t', dtype=str, comments='%', usecols=(0,1)).tolist()
        for i,(name,color) in enumerate(data): # is list of name, color tuples
            # if i>nmax:
            #     break
            mcolors._colors_full_map[name] = color
        custom_colors[category] = data

def register_cmaps():
    """
    Register colormaps and cycles in the cmaps directory.
    Note all of those methods simply modify the dictionary mcm.cmap_d.
    """
    # First read from file
    for file in glob(f'{os.path.dirname(__file__)}/cmaps/*'):
        # Read table of RGB values
        if not re.search('.rgb$', file) and not re.search('.hex$', file):
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
        else:
            cmap = [*open(file)] # just a single line
            if len(cmap)==0:
                continue # file is empty
            cmap = cmap[0].strip().split(',') # csv hex strings
            cmap = np.array([mcolors.to_rgb(c) for c in cmap]) # from list of tuples
        # Register as ListedColormap or LinearSegmentedColormap
        if 'lines' in name.lower():
            cmap   = mcolors.ListedColormap(cmap)
            cmap_r = cmap.reversed() # default name is name+'_r'
            custom_cycles[name] = cmap
        else:
            N = len(cmap) # simple as that; number of rows of colors
            cmap   = mcolors.LinearSegmentedColormap.from_list(name, cmap, N) # using static method is way easier
            cmap_r = cmap.reversed() # default name is name+'_r'
            custom_cmaps[name] = cmap
        # Register maps (this is just what register_cmap does)
        mcm.cmap_d[cmap.name]   = cmap
        mcm.cmap_d[cmap_r.name] = cmap_r
    # Next register names so that they can be invoked ***without capitalization***
    # This always bugged me! Note cannot change dictionary during iteration.
    ignorecase = {}
    for name,cmap in mcm.cmap_d.items():
        if re.search('[A-Z]',name):
            ignorecase[name.lower()] = cmap
    mcm.cmap_d.update(ignorecase)

def register_cycles():
    """
    Register cycles defined right here by dictionaries.
    """
    # Simply register them as ListedColormaps
    for name,colors in cycles.items():
        mcm.cmap_d[name]        = mcolors.ListedColormap([to_rgb(color) for color in colors])
        mcm.cmap_d[f'{name}_r'] = mcolors.ListedColormap([to_rgb(color) for color in colors[::-1]])

# Register stuff when this module is imported
custom_colors = {} # initialize
custom_cycles = {}
custom_cmaps  = {}
register_colors()
register_cmaps()
register_cycles()
print('Registered colors and colormaps.')

#------------------------------------------------------------------------------#
# Generalized colormap/cycle constructors
#------------------------------------------------------------------------------#
def Colormap(*args, light=True, extend='both', ratios=1, resample=False,
        name='custom', register=True, N=None, **kwargs):
    """
    Convenience function for generating colormaps in a variety of ways.
    The 'extend' property will be used to resample LinearSegmentedColormap
    if we don't intend to use both out-of-bounds colors; otherwise we lose
    the strongest colors at either end of the colormap.

    You can still use extend='neither' in Colormap() call with extend='both'
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

    Turns out pcolormesh makes QuadMesh, which itself is a Collection,
    which itself gets colors when calling draw() using update_scalarmappable(),
    which itself uses to_rgba() to get facecolors, which itself is an inherited
    ScalarMappable method that simply calls the colormap with numbers. Anyway
    the issue *has* to be with pcolor, because when giving pcolor an actual
    instance, no longer does that thing where final levels equal extensions.
    Since collection API does nothing to underlying data or cmap, must be
    something done by pcolormesh function.
    """
    N = N or 256 # default
    cmaps = []
    if len(args)==0:
        raise ValueError('Function requires at least 1 positional arg.')
    for cmap in args:
        # Retrieve Colormap instance
        if not isinstance(cmap, mcolors.Colormap):
            if utils.isnumber(cmap):
                cmap = f'C{cmap}' # use current color cycle
            if type(cmap) is dict:
                # Dictionary of hue/sat/luminance values or 2-tuples representing linear transition
                cmap = PerceptuallyUniformColormap.from_hsl(name, **cmap)
                # cmap = space_cmap(name=name, **cmap)
            # (len(cmap)==2 and type(cmap[0]) is str and type(cmap[1]) is dict):
            elif type(cmap) is str:
                # Map name or color for generating monochrome gradiation
                # if type(cmap) is not str:
                #     cmap, cmap_kw = cmap[0], {**cmap_kw, **cmap[1]}
                if cmap in mcm.cmap_d:
                    cmap = mcm.cmap_d[cmap] # get the instance
                else:
                    # Parse extra options
                    cmap_kw = kwargs.copy() # may be different for each cmap in *args
                    regex = '_([rlwdb]+)$'
                    match = re.search(regex, cmap) # declare options with _[flags]
                    cmap = re.sub(regex, '', cmap) # remove options
                    options = '' if not match else match.group(1)
                    if 'r' in options:
                        cmap_kw.update({'reverse':True})
                    if 'w' in options or 'l' in options:
                        light = True
                        if 'w' in options:
                            cmap_kw.update({'white':'white'}) # use *actual* white
                    if 'd' in options or 'b' in options:
                        light = False
                        if 'b' in options:
                            cmap_kw.update({'black':'black'}) # use *actual* black
                    # Build colormap
                    cmap = to_rgb(cmap) # to ensure is hex code/registered color
                    if light:
                        cmap = light_cmap(cmap, name=name, **cmap_kw)
                    else:
                        cmap = dark_cmap(cmap, name=name, **cmap_kw)
            else:
                # List of colors
                cmap = mcolors.ListedColormap(cmap, name=name)
        cmaps += [cmap]
    cmap = merge_cmaps(cmaps, name=name, ratios=ratios, **kwargs)
    name = cmap.name # just to make sure
    # Optionally make extremes same color as map
    offset = {'neither':-1, 'max':0, 'min':0, 'both':1}
    if extend not in offset:
        raise ValueError(f'Unknown extend option {extend}.')
    if isinstance(cmap, mcolors.LinearSegmentedColormap):
        if resample:
            cmap = cmap._resample(N-offset[extend]) # see mcm.get_cmap source
        else:
            pass
            # cmap = pad_cmap(cmap, 1/float(N), extend=extend)
    # Optionally register a colormap
    if name and register:
        if name.lower() in [cat_cmap.lower() for cat,cat_cmaps in categories_default.items()
                    for cat_cmap in cat_cmaps if 'Custom' not in cat]:
            raise ValueError(f'Builtin colormap "{name}" already exists. Choose a different name.')
        elif name in mcm.cmap_d:
            pass # no warning necessary
            # print(f'Warning: Overwriting existing colormap "{name}".')
        mcm.cmap_d[name] = cmap
        if re.search('[A-Z]',name):
            mcm.cmap_d[name.lower()] = cmap
        # print(f'Registered name {name}.') # not necessary
    return cmap

def Cycle(*args, vmin=0, vmax=1):
    """
    Convenience function to draw colors from arbitrary ListedColormap or
    LinearSegmentedColormap. Use vmin/vmax to scale your samples.
    """
    if len(args)==0:
        raise ValueError('Function requires at least 1 positional arg.')
    if len(args)==1:
        samples = 10
    else:
        args, samples = args[:-1], args[-1]
    cmap = Colormap(*args) # the cmap object itself
    if isinstance(cmap, mcolors.ListedColormap):
        # Just get the colors
        colors = cmap.colors
    elif isinstance(cmap, mcolors.LinearSegmentedColormap):
        # Employ ***more flexible*** version of get_cmap() method, which does this:
        # LinearSegmentedColormap(self.name, self._segmentdata, lutsize)
        if utils.isnumber(samples):
            # samples = np.linspace(0, 1-1/nsample, nsample) # from 'centers'
            samples = np.linspace(0, 1, samples) # from edge to edge
        else:
            samples = np.array(samples)
        colors = cmap((samples-vmin)/(vmax-vmin))
    else:
        raise ValueError(f'Colormap returned weird object type: {type(cmap)}.')
    return colors

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
    def __init__(self, name, segmentdata, space='hsl', scale=True, mask=False, **kwargs):
        """
        Initialize with dictionary of values.
        Arguments
        ---------
            scale : The input hues, saturations, and luminances should all be
                normalized to the range 0-1. To disable this behavior, set
                scale=False.
            mask : Whether to mask out-of-range colors as black, or just clip
                the RGB values (distinct from colormap clipping the extremes).
        """
        # Attributes
        space = get_space(space)
        self.space = space
        self.mask  = mask
        # First sanitize the segmentdata by converting color strings to their
        # corresponding channel values
        keys   = {*segmentdata.keys()}
        target = {'hue', 'saturation', 'luminance'}
        if keys != target and keys != {*target, 'alpha'}:
            raise ValueError('Invalid segmentdata dictionary.')
        for key,array in segmentdata.items():
            for i,xyy in enumerate(array):
                xyy = list(xyy) # make copy!
                for j,y in enumerate(xyy[1:]):
                    j += 1 # fix
                    xyy[j] = get_channel(y, key, space, scale) # also *scales* values to 99, 99, 359
                segmentdata[key][i] = xyy
        # Initialize
        super().__init__(name, segmentdata, **kwargs)
        # self.cdict = dict(hue={}, saturation={}, luminance={}) # to minimize string-name lookup

    def _init(self):
        """
        As with LinearSegmentedColormap, but convert each value
        in the lookup table from 'input' to RGB.
        """
        # First generate the lookup table
        scale = get_scale(self.space)
        self._lut = np.ones((self.N+3, 4), float) # fill
        for i,key in enumerate(('hue','saturation','luminance')):
            array       = np.array(self._segmentdata[key])
            array[:,1:] = array[:,1:]
            # /scale[i] # scale the y-values
            # if key=='hue':
            #     print('initial hue array', array)
            #     array[:,1:] = array[:,1:] % 1 # allow circular hues
            self._lut[:-3,i] = make_mapping_array(self.N, array, self._gamma)
        if 'alpha' in self._segmentdata:
            self._lut[:-3,3] = make_mapping_array(self.N, self._segmentdata['alpha'], 1)
        # Make hues circular
        self._lut[:-3,0] %= 359 # mod
        # if key=='hue':
        #     print('initial hue array', array)
        #     array[:,1:] = array[:,1:] % 1 # allow circular hues
        # Set extremes
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
        return PerceptuallyUniformColormap(self.name, self._segmentdata, scale=False, space=self.space, N=N)

    @staticmethod
    def from_hsl(name, h=1.0, s=1.0, l=[1,0.2], c=None, a=None, **kwargs):
        """
        Simply wrapper for from_list. Constructs list of colors from user
        input h=hues, s=saturations, l=luminances.
        """
        # Get channels from keyword arguments
        if c is not None:
            s = c
        if a is None:
            a = 1.0
        nlevs = {len(c) for c in (s,h,l) if utils.isvector(c)} # set
        if 1 in nlevs:
            nlevs.remove(1)
        if len(nlevs)==0:
            raise ValueError('At least one channel level specifier has to be non-scalar.')
        elif len(nlevs)!=1:
            raise ValueError(f'Got {len(h)} hue, {len(s)} saturation, {len(l)} luminance values.')
        nlevs = nlevs.pop()
        nsegs = nlevs-1
        hues   = [h]*nlevs if utils.isscalar(h) else [*h]
        sats   = [s]*nlevs if utils.isscalar(s) else [*s]
        lums   = [l]*nlevs if utils.isscalar(l) else [*l]
        alphas = [a]*nlevs if utils.isscalar(a) else [*a]
        # Now build hsl tuples and pass
        colors = [[h,s,l,a] for h,s,l,a in zip(hues,sats,lums,alphas)]
        return PerceptuallyUniformColormap.from_list(name, colors, **kwargs)

    @staticmethod
    def from_list(name, colors, ratios=None, **kwargs):
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
        # Check input
        if not np.iterable(colors):
            raise TypeError('Colors must be iterable.')
        if ratios is not None:
            xvals = np.atleast_1d(ratios) # could be ratios=1, i.e. dummy
            if len(xvals) != len(colors) - 1:
                raise ValueError(f'Got {len(colors)} colors, but {len(ratios)} ratios.')
            xvals = np.concatenate(([0], np.cumsum(xvals)))
            xvals = xvals/np.max(xvals) # normalize to 0-1
        else:
            xvals = np.linspace(0,1,len(colors))
        # Build dictionary
        cdict = dict(hue=[], saturation=[], luminance=[], alpha=[])
        for x,color in zip(xvals,colors):
            # Assign after ensuring alpha channel exists
            color = add_alpha(color)
            h, s, l, a = color # get values
            cdict['hue'].append((x, h, h))
            cdict['saturation'].append((x, s, s))
            cdict['luminance'].append((x, l, l))
            cdict['alpha'].append((x, a, a))
        return PerceptuallyUniformColormap(name, cdict, **kwargs)

def make_mapping_array(N, data, gamma=1.0):
    """
    Carbon copy of matplotlib version, but this one doesn't clip values.
    """
    # Get array
    try:
        adata = np.array(data)
    except Exception:
        raise TypeError("data must be convertible to an array")
    shape = adata.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError("data must be nx3 format")
    # Get indices
    x  = adata[:, 0]
    y0 = adata[:, 1]
    y1 = adata[:, 2]
    if x[0] != 0.0 or x[-1] != 1.0:
        raise ValueError('Data mapping points must start with x=0 and end with x=1')
    if (np.diff(x) < 0).any():
        raise ValueError('Data mapping points must have x in increasing order')
    # Begin generation of lookup table
    x    = x * (N - 1)
    lut  = np.zeros((N,), float)
    xind = (N - 1) * np.linspace(0, 1, N) ** gamma
    ind  = np.searchsorted(x, xind)[1:-1]
    # Calculate
    distance = (xind[1:-1] - x[ind - 1]) / (x[ind] - x[ind - 1])
    lut[1:-1] = distance * (y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0] = y1[0]
    lut[-1] = y0[-1]
    return lut

#------------------------------------------------------------------------------#
# Colormap constructors
#------------------------------------------------------------------------------#
def merge_cmaps(cmaps, n=512, name='merged', ratios=1, **kwargs):
    """
    Merge arbitrary colormaps.
    Arguments
    ---------
        cmaps : 
            List of colormap strings or instances for merging.
        n :
            Number of lookup table colors desired for output colormap.
        name :
            Name of output colormap.
    Notes
    -----
    Old method had us simply calling the colormap with arrays of fractions.
    This was sloppy, because it just samples locations on the lookup table and
    will therefore degrade the original, smooth, functional transitions.
    Better method is to combine the _segmentdata arrays and simply scale
    the x coordinates in each (x,y1,y2) channel-tuple according to the ratios.
    In the case of ListedColormaps, we just combine the colors.
    """
    # Return
    if len(cmaps)==1:
        return cmaps[0] # needed to avoid recursion!
    # Parse ratios
    ratios = ratios or 1
    if utils.isscalar(ratios):
        ratios = [1]*len(cmaps)
    # Merge generalized colormap objects
    # colors = [color for cmap,ratio in zip(cmaps,ratios) for color in
    #           Colormap(cmap, **kwargs)(np.linspace(0,1,int(n*ratio)))]
    # return smooth_cmap(colors, n=n, name=name)
    # More accurate method is below
    cmaps = [Colormap(cmap, N=None, **kwargs) for cmap in cmaps] # set N=None to disable resamping
    if all(isinstance(cmap,mcolors.ListedColormap) for cmap in cmaps):
        if not np.all(ratios==1):
            raise ValueError(f'Cannot assign different ratios when mering ListedColormaps.')
        colors = [color for cmap.colors in cmaps for color in cmap.colors]
        cmap = mcolors.ListedColormap(colors, name=name, N=len(colors))
    elif all(isinstance(cmap,mcolors.LinearSegmentedColormap) for cmap in cmaps):
        kinds = {type(cmap) for cmap in cmaps}
        if len(kinds)>1:
            raise ValueError(f'Got mixed colormap types.')
        kind = kinds.pop() # colormap kind
        keys = {key for cmap in cmaps for key in cmap._segmentdata.keys()}
        if len(keys) not in (3,4):
            raise ValueError(f'Got mixed segmentdata keys: {keys}.')
        ratios = np.array(ratios)/np.sum(ratios) # so if 4 cmaps, will be 1/4
        coords = np.concatenate([[0], np.cumsum(ratios)])
        widths = coords[1:] - coords[:-1]
        segmentdata = {}
        for key in keys:
            datas = []
            for x,width,cmap in zip(coords[:-1],widths,cmaps):
                data = np.array(cmap._segmentdata[key])
                data[:,0] = x + width*data[:,0]
                datas.append(data)
            data = np.concatenate(datas, axis=0)
            data[:,0] = data[:,0]/data[:,0].max(axis=0) # scale to make maximum 1
            segmentdata[key] = data
        kwargs = dict(N=n)
        if kind is PerceptuallyUniformColormap:
            spaces = {cmap.space for cmap in cmaps}
            if len(spaces)>1:
                raise ValueError(f'Trying to merge colormaps with different HSL spaces {repr(spaces)}.')
            kwargs.update(dict(space=spaces.pop(), scale=False))
        cmap = kind(name, segmentdata, **kwargs)
    else:
        raise ValueError('All colormaps should be of the same type (Listed or LinearSegmented).')
    return cmap

def light_cmap(color, reverse=False, space='hsl', white='#eeeeee', light=None, **kwargs):
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
    white = light or white
    space = get_space(space)
    _, ws, wl = to_xyz(to_rgb(white), space)
    h, s, l = to_xyz(to_rgb(color), space)
    h, s, l, ws, wl = h/359.0, s/99.0, l/99.0, ws/99.0, wl/99.0
    ws = s # perhaps don't change the chroma by default
    index = slice(None,None,-1) if reverse else slice(None)
    return PerceptuallyUniformColormap.from_hsl(h, [s,ws][index], [l,wl][index], space=space, **kwargs)
    # return space_cmap(n, h, [s,ws][index], [l,wl][index], space=space, **kwargs)

def dark_cmap(color, reverse=False, space='hsl', black='#444444', dark=None, gray=None, grey=None, **kwargs):
    """
    Make a sequential colormap that blends from gray to color.
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
    black = grey or gray or dark or black # alternative kwargs
    space = get_space(space)
    _, bs, bl = to_xyz(to_rgb(black), space)
    h, s, l = to_xyz(to_rgb(color), space)
    h, s, l, bs, bl = h/359.0, s/99.0, l/99.0, bs/99.0, bl/99.0
    index = slice(None,None,-1) if reverse else slice(None)
    return PerceptuallyUniformColormap.from_hsl(h, [bs,s][index], [bl,l][index], space=space, **kwargs)
    # return space_cmap(n, h, [bs,s][index], [bl,l][index], space=space, **kwargs)

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
                    print(f'Warning: {message} channel {n} (<0).')
            if c>1:
                if mask:
                    rgb[-1] = 1 # full alpha
                    for j in range(3):
                        rgb[j] = gray
                else:
                    rgb[i] = 1 # keep
                if n not in over:
                    over[n] = True
                    print(f'Warning: {message} channel {n} (>1).')
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
    colors = Cycle(cmap, samples)
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
    if cycle=='reset':
        colors = [(0.0, 0.0, 1.0), (0.0, .50, 0.0), (1.0, 0.0, 0.0), (.75, .75, 0.0),
                  (.75, .75, 0.0), (0.0, .75, .75), (0.0, 0.0, 0.0)]
    elif cycle in seaborn_cycles:
        colors = cycles[cycle] + [(.1, .1, .1)]
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
        norm = 'step' # the default case when levels are provided
    elif norm is None:
        norm = mcolors.Normalize() # default is just linear from 0 to 1
    elif type(norm) is not str: # dictionary lookup
        raise ValueError(f'Unknown norm "{norm}".')
    if type(norm) is str:
        if norm not in normalizers:
            raise ValueError(f'Unknown normalizer "{norm}". Options are {", ".join(normalizers.keys())}.')
        norm = normalizers[norm](**kwargs)
    return norm

class LinearSegmentedNorm(mcolors.Normalize):
    """
    As in BoundaryNorm case, but instead we linearly *interpolate* colors
    between the provided boundary indices. Use this e.g. if you want to
    have the gradation from 0-1 'zoomed in' and gradation from 0-10 'zoomed out'.
    In this case, can control number of colors by setting the *lookup table*
    when declaring colormap.
    """
    def __init__(self, levels, midpoint=None, clip=False, ncolors=None, **kwargs):
        # Very simple
        levels = np.atleast_1d(levels)
        if np.any((levels[1:]-levels[:-1])<=0):
            raise ValueError(f'Levels passed to LinearSegmentedNorm must be monotonically increasing.')
        super().__init__(np.nanmin(levels), np.nanmax(levels), clip) # second level superclass
        self.levels = levels # alias for boundaries

    def __call__(self, value, clip=None):
        # TODO: Add optional midpoint, this class will probably end up being one of
        # my most used if so. Midpoint would just ensure <value> corresponds to 0.5 in cmap
        # Map data values in range (vmin,vmax) to color indices in colormap
        # If have 11 levels and between ids 9 and 10, interpolant is between .9 and 1
        value = np.atleast_1d(value)
        norm  = np.empty(value.shape)
        for i,v in enumerate(value.flat):
            if np.isnan(v):
                continue
            idxs = np.linspace(0, 1, self.levels.size)
            locs = np.where(v>=self.levels)[0]
            if locs.size==0:
                norm[i] = 0
            elif locs.size==self.levels.size:
                norm[i] = 1
            else:
                cutoff      = locs[-1]
                interpolant = idxs[cutoff:cutoff+2] # the 0-1 normalization
                interpolee  = self.levels[cutoff:cutoff+2] # the boundary level values
                norm[i]     = np.interp(v, interpolee, interpolant) # interpolate between 2 points
        return ma.masked_array(norm, np.isnan(value))

    def inverse(self, norm):
        # Performs inverse operation of __call__
        norm  = np.atleast_1d(norm)
        value = np.empty(norm.shape)
        for i,n in enumerate(norm.flat):
            if np.isnan(n):
                continue
            idxs = np.linspace(0, 1, self.levels.size)
            locs = np.where(n>=idxs)[0]
            if locs.size==0:
                value[i] = np.nanmin(self.levels)
            elif locs.size==self.levels.size:
                value[i] = np.nanmax(self.levels)
            else:
                cutoff = locs[-1]
                interpolee  = idxs[cutoff:cutoff+2] # the 0-1 normalization
                interpolant = self.levels[cutoff:cutoff+2] # the boundary level values
                value[i]    = np.interp(n, interpolee, interpolant) # interpolate between 2 points
        return ma.masked_array(value, np.isnan(norm))

class StepNorm(mcolors.BoundaryNorm):
    """
    Simple normalizer that *interpolates* from an RGB array at point
    (level_idx/num_levels) along the array, instead of choosing color
    from (transform(level_value)-transform(vmin))/(transform(vmax)-transform(vmin))
    where transform can be linear, logarithmic, etc.

    Note
    ----
    If you are using a diverging colormap with extend='max/min', the center
    will get messed up. But that is very strange usage anyway... so you'll
    probably be fine.

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
    def __init__(self, levels, centers=False, clip=False, ncolors=None, extend=None, **kwargs):
        # Don't need to call partent initializer, this is own implementation; just need
        # it to be subclass so ColorbarBase methods will detect it
        # See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        extend = extend or 'both'
        levels = np.atleast_1d(levels)
        if np.any((levels[1:]-levels[:-1])<=0):
            raise ValueError('Levels passed to Normalize() must be monotonically increasing.')
        if extend not in ('both','min','max','neither'):
            raise ValueError(f'Unknown extend option "{extend}". Choose from "min", "max", "both", "neither".')
        # Bins were passed, not edges
        if centers:
            levels = utils.edges(levels)
        # Custom attributes
        offset = {'both':2, 'min':1, 'max':1, 'neither':0}
        slices = {'both':slice(1,-1),   'min':slice(1,None),
                  'max':slice(None,-1), 'neither':slice(None)}
        self.levels   = levels
        self.normvals = np.linspace(0, 1, offset[extend]+levels.size-1)[slices[extend]]
        # Required
        self.N = len(levels)
        self.boundaries = levels # alias read by other functions
        self.vmin = levels[0]
        self.vmax = levels[-1]
        self.clip = clip

    def __call__(self, value, clip=None):
        # TODO: Add optional midpoint, this class will probably end up being one of
        # my most used if so. Midpoint would just ensure <value> corresponds to 0.5 in cmap
        # Map data values in range (vmin,vmax) to color indices in colormap
        # If have 11 levels and between ids 9 and 10, interpolant is between .9 and 1
        value = np.atleast_1d(value)
        norm  = np.empty(value.shape)
        for i,v in enumerate(value.flat):
            if np.isnan(v):
                continue
            locs = np.where(v>=self.levels)[0]
            if locs.size==0:
                norm[i] = 0
            elif locs.size==self.levels.size:
                norm[i] = 1
            else:
                cutoff  = locs[-1]
                norm[i] = self.normvals[cutoff]
                # (self.idxs[cutoff] + self.idxs[cutoff+1])/2
        return ma.masked_array(norm, np.isnan(value))

    def inverse(self, norm):
        # Not possible
        raise ValueError('StepNorm is not invertible.')

class StretchNorm(mcolors.Normalize):
    """
    Class that can 'stretch' and 'compress' either side of a colormap about
    some midpoint, proceeding exponentially (exp>0) or logarithmically (exp<0)
    down the linear colormap from the center point. Default midpoint is vmin, i.e.
    we just stretch to the right. For diverging colormaps, use midpoint 0.5.
    """
    def __init__(self, exp=0, extend='neither', midpoint=None, vmin=None, vmax=None, clip=None):
        # User will use -10 to 10 scale; converted to value used in equation
        if abs(exp) > 10: raise ValueError('Warping scale must be between -10 and 10.')
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        self.exp = exp
        self.extend = extend
        # mcolors.Normalize.__init__(self, vmin, vmax, clip)

    # Function
    def warp(x, exp, exp_max=4):
        # Returns indices stretched so neutral/low values are sampled more heavily
        # Will artifically use exp to signify stretching away from neutral vals,
        # or compressing toward neutral vals
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
        # Initial stuff
        if self.midpoint is None:
            midpoint = self.vmin
        else:
            midpoint = self.midpoint
        # Get middle point in 0-1 coords, and value
        midpoint_scaled = (midpoint - self.vmin)/(self.vmax - self.vmin)
        value_scaled    = (value - self.vmin)/(self.vmax - self.vmin)
        try: iter(value_scaled)
        except TypeError:
            value_scaled = np.arange(value_scaled)
        value_cmap = ma.empty(value_scaled.size)
        for i,v in enumerate(value_scaled):
            # Get values, accounting for midpoints
            if v<0:
                v = 0
            if v>1:
                v = 1
            if v>=midpoint_scaled:
                block_width = 1 - midpoint_scaled
                value_cmap[i] = (midpoint_scaled + 
                    block_width*self.warp((v - midpoint_scaled)/block_width, self.exp)
                    )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled - 
                        block_width*self.warp((midpoint_scaled - v)/block_width, self.exp)
                        )
        if self.extend=='both' or self.extend=='max':
            value_cmap[value_cmap>1] = 1
        if self.extend=='both' or self.extend=='min':
            value_cmap[value_cmap<0] = 0
        return value_cmap

#------------------------------------------------------------------------------#
# Visualizations
#------------------------------------------------------------------------------#
def color_show(groups=['open',['crayons','xkcd']], ncols=4, nbreak=15, minsat=0.1):
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
        group = group or 'open'
        if type(group) is str:
            group = [group]
        colors = {}
        for name in group:
            # Read colors from current cycler
            if name=='cycle':
                seen = set() # trickery
                cycle_colors = rcParams['axes.prop_cycle'].by_key()['color']
                cycle_colors = [color for color in cycle_colors if not (color in seen or seen.add(color))] # trickery
                colors.update({f'C{i}':v for i,v in enumerate(cycle_colors)})
            # Read custom defined colors
            else:
                colors.update({name:color for name,color in custom_colors[name]}) # convert from list of length-2 tuples to dictionary
        # Group colors together by discrete range of hue, then sort by value
        # For opencolors this is not necessary
        if 'open' in group:
            space = 0.5
            swatch = 1.5
            names = ['gray', 'red', 'pink', 'grape', 'violet', 'indigo', 'blue', 'cyan', 'teal', 'green', 'lime', 'yellow', 'orange']
            nrows, ncols = 10, len(names) # rows and columns
            sorted_names = [[name+str(i) for i in range(nrows)] for name in names]
        # For other palettes this is necessary
        else:
            # Keep in separate columns
            space = 1
            swatch = 1
            ncols = nbreak-1 # group by breakpoint
            colors_hsv = {k:tuple(colormath.rgb_to_hsv(*mcolors.to_rgb(v))) for k,v in colors.items()}
            breakpoints = np.linspace(0,1,nbreak) # group in blocks of 20 hues
            sorted_names = [] # initialize
            testsat = (lambda x: x<minsat) # test saturation
            for n in range(len(breakpoints)):
                if n==0: # grays
                    fcolors = [(name,hsv) for name,hsv in colors_hsv.items()
                        if testsat(hsv[1])]
                    sortfunc = lambda x: 3*x[2]+x[0]
                else: # colors
                    start, end = breakpoints[n-1], breakpoints[n]
                    testhue = (lambda x: start<=x<=end) if end is breakpoints[-1] \
                        else (lambda x: start<=x<end) # two possible tests
                    fcolors = [(name,hsv) for name,hsv in colors_hsv.items()
                        if testhue(hsv[0]) and not testsat(hsv[1])] # grays have separate category
                    sortfunc = lambda x: 3*x[2]+x[1]
                sorted_index = np.argsort([sortfunc(v[1]) for v in fcolors]) # indices to build sorted list
                sorted_names.append([fcolors[i][0] for i in sorted_index]) # append sorted list
            nrows = max(len(huelist) for huelist in sorted_names) # number of rows
            # Now concatenate those columns so get nice rectangle
            names = [i for sublist in sorted_names for i in sublist]
            sorted_names = [[]]
            nrows = len(names)//ncols+1
            for i,name in enumerate(names):
                if ((i + 1) % nrows)==0:
                    sorted_names.append([]) # add new empty list
                sorted_names[-1].append(name)
        # Create plot by iterating over columns 
        figsize = (8*space*(ncols/4), 5*(nrows/40)) # 5in tall with 40 colors in column
        fig, ax = plt.subplots(figsize=figsize)
        X, Y = fig.get_dpi()*fig.get_size_inches() # size in *dots*; make these axes units
        h, w = Y/(nrows+1), X/ncols # height and width of row/column in *dots*
        for col,huelist in enumerate(sorted_names):
            for row,name in enumerate(huelist): # list of colors in hue category
                y = Y - (row * h) - h
                xi_line = w*(col + 0.05)
                xf_line = w*(col + 0.25*swatch)
                xi_text = w*(col + 0.25*swatch + 0.03*swatch)
                print_name = name.split('xkcd:')[-1] # make sure no xkcd:
                ax.text(xi_text, y, print_name, fontsize=h*0.8, ha='left', va='center')
                ax.hlines(y+h*0.1, xi_line, xf_line, color=colors[name], lw=h*0.6)
        ax.set_xlim(0,X)
        ax.set_ylim(0,Y)
        ax.set_axis_off()
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)
        # Save
        fig.savefig(f'{os.path.dirname(__file__)}/colors_{"-".join(group)}.pdf',
                bbox_inches='tight', format='pdf')
        figs += [fig]
    return figs

def cycle_show():
    """
    Show off the different color cycles.
    Wrote this one myself, so it uses the custom API.
    """
    # cycles = plt.get_cycles() # function should have been added by the rc plugin
    cycles = mcolors.CYCLES
    nrows = len(cycles)//2+len(cycles)%2
    fig, axs = plt.subplots(figsize=(6,nrows*1.5), ncols=2, nrows=nrows)
    axs = [ax for sub in axs for ax in sub]
    fig.subplots_adjust(top=.99, bottom=.01, left=.02, right=0.98, hspace=.2, wspace=.02)
    state = np.random.RandomState(123412)
    for i,(ax,(key,value)) in enumerate(zip(axs,cycles.items())):
        propcycle = cycler('color', value)
        ax.set_prop_cycle(propcycle)
        lines = ax.plot(state.rand(10,len(value)), lw=5, ls='-')
        for j,l in enumerate(lines):
            l.set_zorder(len(lines)-j) # make first lines have big zorder
        title = f'{key}: {len(value)} colors'
        ax.set_xlim((-0.5,10))
        ax.set_title(title)
        for axis in 'xy':
            ax.tick_params(axis=axis, which='both', labelbottom=False, labelleft=False,
                    bottom=False, top=False, left=False, right=False)
    if len(cycles)%2==1:
        axs[-1].set_visible(False)
    # Save
    fig.savefig(f'{os.path.dirname(__file__)}/cycles.pdf',
            bbox_inches='tight', format='pdf')
    return fig

# def cmap_show(N=31, ignore=['Miscellaneous','Sequential2','Diverging2']):
def cmap_show(N=31, ignore=['Diverging Alt', 'Sequential Alt', 'Rainbow Alt', 'Miscellaneous']):
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
    cmaps_all = [cmap for cmap in mcm.cmap_d.keys() if
            not cmap.endswith('_r')
            and not re.search('[A-Z]',cmap)
            and 'Vega' not in cmap
            and (isinstance(mcm.cmap_d[cmap],mcolors.LinearSegmentedColormap) or cmap in exceptions)]
    cmaps_listed = [cmap for cmap in mcm.cmap_d.keys() if
            not cmap.endswith('_r')
            and not re.search('[A-Z]',cmap)
            and 'Vega' not in cmap
            and (not isinstance(mcm.cmap_d[cmap],mcolors.LinearSegmentedColormap) and cmap not in exceptions)]
    # Detect unknown/manually created colormaps, and filter out
    # colormaps belonging to certain section
    cmaps_ignore  = [cmap for cat,cmaps in categories_default.items() for cmap in cmaps if cat in ignore]
    categories    = {cat:cmaps for cat,cmaps in categories_default.items() if cat not in ignore}
    cmaps_known   = [cmap for cat,cmaps in categories.items() for cmap in cmaps if cmap in cmaps_all]
    cmaps_missing = [cmap for cat,cmaps in categories.items() for cmap in cmaps if cmap not in cmaps_all]
    cmaps_custom  = [cmap for cmap in cmaps_all if cmap not in cmaps_known and cmap not in cmaps_ignore]
    # Attempt to auto-detect diverging colormaps, just sample the points on either end
    # Do this by simply summing the RGB channels to get HSV brightness
    categories['Custom Diverging'][:] = [] # empty list out again
    categories['Custom Sequential'][:] = []
    for cmap in cmaps_custom:
        m = mcm.cmap_d[cmap]
        l = lambda i: to_xyz(to_rgb(m(i)), 'hcl')[2] # get luminance
        if (l(0)<l(0.5) and l(1)<l(0.5)): # or (l(0)>l(0.5) and l(1)>l(0.5)):
            categories['Custom Diverging'] += [cmap]
        else:
            categories['Custom Sequential'] += [cmap]
    # categories['Custom'] = cmaps_custom
    # categories['Custom'] = [cmap for cmap in cmaps_custom if cmap[:2]!='cb']
    # categories['ColorBrewer'] = [cmap for cmap in cmaps_custom if cmap[:2]=='cb']
    if cmaps_missing:
        print(f'Missing colormaps: {", ".join(cmaps_missing)}')
    if cmaps_ignore:
        print(f'Ignored colormaps: {", ".join(cmaps_ignore)}')
    # print(f'Custom colormaps: {", ".join(cmaps_custom)}')
    # print(f'Listed colormaps: {", ".join(cmaps_listed)}')
    # Attempt sorting based on hue
    for cat in ['Custom Sequential', 'Builtin Sequential']:
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
    fig.subplots_adjust(top=0.98, bottom=0.01, left=0.15, right=0.99)
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
    filename = f'{os.path.dirname(__file__)}/colormaps.pdf'
    print(f"Saving figure to: {filename}.")
    fig.savefig(filename, bbox_inches='tight')
    return fig

# Finally our dictionary of normalizers
# Includes some custom classes, so has to go at end
normalizers = {
    'none':       mcolors.NoNorm,
    'null':       mcolors.NoNorm,
    'step':       StepNorm,
    'segmented':  LinearSegmentedNorm,
    'boundary':   mcolors.BoundaryNorm,
    'log':        mcolors.LogNorm,
    'linear':     mcolors.Normalize,
    'power':      mcolors.PowerNorm,
    'symlog':     mcolors.SymLogNorm,
    }

# TODO: Figure out if this is still useful
# ***Inverse of cycle_factory.***
# Generate colormap instance from list of levels and colors.
# * Generally don't want these kinds of colorbars to 'extend', but you can.
# * Object will assume one color between individual levels; for example,
#   if levels are [0,0.5,0.6,1], the first and last color-zones will be way thicker.
# * If you want unevenly space intervals but evenly spaced boundaries, use custom
#   Norm instead of default.
# * This is *one way* to create a colormap from list of colors; another way is
#   to just pass a list of colors to any _cformat method. That method will
#   also generate levels that are always equally spaced.
# if levels is None:
#     levels = np.linspace(0,1,len(colors)+1)
# if len(levels)!=len(colors)+1:
#     raise ValueError(f"Have {len(levels):d} levels and {len(colors):d} colors. Need ncolors+1==nlevels.")
# if extend=='min':
#     colors = ['w', *colors] # add dummy color
# elif extend=='max':
#     colors = [*colors, 'w'] # add dummy color
# elif extend=='both':
#     colors = ['w', *colors, 'w']
# elif extend!='neither':
#     raise ValueError("Unknown extend option \"{extend}\".")
# cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend=extend) # creates ListedColormap!!!
# return cmap
# Generate mappable using contourf
# def Mappable(cmap, samples=None, levels=None, norm=None):
# cmap = cmap_factory(cmap)
# if norm is None:
#     norm = mcolors.Normalize(0,1)
# elif not hasattr(norm,'__call__'):
#     raise ValueError('Norm has to be callable.')
# elif not isinstance(norm,mcolors.Normalize):
#     raise ValueError('Norm has to be an mcolors.Normalize class or subclass.')
# # Specify the *centers*, then interpolate to *edges* according to normalizer
# if samples is not None:
#     samples = np.array(samples)
#     if not isinstance(cmap, mcolors.BoundaryNorm):
#         samples = norm(samples) # translate to normalized-space (e.g. LogNorm)
#     if samples[1]<samples[0]:
#         samples = samples[::-1]
#         reverse = True
#     idelta = samples[1]-samples[0]
#     fdelta = samples[-1]-samples[-2]
#     levels = [samples[0]-idelta/2, *((samples[1:]+samples[:-1])/2), samples[-1]+fdelta/2]
#     if not isinstance(cmap, mcolors.BoundaryNorm):
#         levels = norm.inverse(levels) # translate back from normalized-space
# # Optionally reverse the colormap
# if levels[1]<levels[0]: # specify the *levels*
#     levels = levels[::-1]
#     reverse = True
# if reverse:
#     cmap = cmap.reversed() # reverse the cmap
# m = plt.contourf([0,0], [0,0], np.nan*np.ones((2,2)), cmap=cmap, levels=levels, norm=norm)
# return m
# This was dumb, just needed to edit normalizer
# def pad_cmap(cmap, amount, extend='both'):
#     """
#     Pads the extremes for lookup table in cmap, so that, if we do not ever want
#     to color our 'extremes', we can still draw levels that sample the full range
#     of possible colors.
#     The 'amount' is fraction relative to the entire 0-1 width.
#     """
#     # Check stuff
#     if extend not in ('min','max','neither'):
#         raise ValueError(f'Extend must be "min", "max", or "neither".')
#     if not isinstance(cmap, mcolors.LinearSegmentedColormap):
#         raise TypeError('Input must be LinearSegmentedColormap.')
#     # Apply padding
#     offset = 0
#     name = cmap.name
#     segmentdata = cmap._segmentdata
#     if extend in ('max','neither'):
#         offset += amount # pad bottom of colormap
#     for channel in ('red','green','blue'):
#         data = np.array(segmentdata[channel])
#         if extend in ('max','neither'): # pad low values
#         data = np.concatenate((data, pad), axis=1)
#     return mcolors.LinearSegmentedColormap(name, segmentdata, N=segmentdata.shape[0])

