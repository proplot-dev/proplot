#!/usr/bin/env python3
"""
Registers colormaps, color cycles, and color string names with `register_cmaps`,
`register_cycles`, and `register_colors`. Defines the `Colormap` and `Cycle`
tools for creating new colormaps and color cycles. Defines helpful new
`~matplotlib.colors.Normalize` and `~matplotlib.colors.Colormap` classes.
Adds tools for visualizing colorspaces, colormaps, color names, and color
cycles.

See the :ref:`Color usage` section of "Getting Started" for details.
"""
# Potential bottleneck, loading all this stuff?  *No*. Try using @timer on
# register functions, turns out worst is colormap one at 0.1 seconds. Just happens
# to be a big package, takes a bit to compile to bytecode then import.
import os
import re
import json
import glob
import cycler
from lxml import etree
from numbers import Number, Integral
import numpy as np
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib as mpl
import warnings
from . import colormath
from .utils import _notNone, _timer
rcParams = mpl.rcParams
__all__ = [
    'BinNorm', 'CmapDict', 'ColorCacheDict',
    'LinearSegmentedNorm', 'MidpointNorm', 'PerceptuallyUniformColormap',
    'Colormap', 'Cycle', 'Norm',
    'cmaps', 'cycles', 'colordict',
    'fonts', 'fonts_system', 'fonts_proplot',
    'colors',
    'make_mapping_array', 'monochrome_cmap',
    'register_cmaps', 'register_colors', 'register_cycles', 'register_fonts',
    'saturate', 'shade', 'show_cmaps', 'show_channels',
    'show_colors', 'show_colorspaces', 'show_cycles', 'show_fonts',
    'to_rgb', 'to_xyz',
    ]

# Data diretories
DATA_USER = os.path.join(os.path.expanduser('~'), '.proplot')
DATA_USER_CMAPS = os.path.join(DATA_USER, 'cmaps')
DATA_USER_CYCLES = os.path.join(DATA_USER, 'cycles')
DATA_USER_FONTS = os.path.join(DATA_USER, 'fonts') # user fonts
DATA_CMAPS = os.path.join(os.path.dirname(__file__), 'cmaps') # or parent, but that makes pip install distribution hard
DATA_CYCLES = os.path.join(os.path.dirname(__file__), 'cycles') # or parent, but that makes pip install distribution hard
DATA_COLORS = os.path.join(os.path.dirname(__file__), 'colors') # or parent, but that makes pip install distribution hard
DATA_FONTS = os.path.join(os.path.dirname(__file__), 'fonts') # proplot fonts
if not os.path.isdir(DATA_USER):
    os.mkdir(DATA_USER)
if not os.path.isdir(DATA_USER_CMAPS):
    os.mkdir(DATA_USER_CMAPS)
if not os.path.isdir(DATA_USER_CYCLES):
    os.mkdir(DATA_USER_CYCLES)
if not os.path.isdir(DATA_USER_FONTS):
    os.mkdir(DATA_USER_FONTS)

# Colormap stuff
CMAPS_CATEGORIES = {
    # User cmaps, placed here so it always appears on top
    'User': (),
    # Assorted origin, but these belong together
    'Grayscale': (
        'Grays', 'Mono', 'GrayCycle',
        ),
    # Builtin
    'Matplotlib Originals': (
        'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'twilight', 'twilight_shifted',
        ),
    # seaborn
    'Seaborn Originals': (
        'Rocket', 'Mako', 'IceFire', 'Vlag',
        ),
    # PerceptuallyUniformColormap
    'ProPlot Sequential': (
        'Fire',
        'Stellar',
        'Boreal',
        'Marine',
        'Dusk',
        'Glacial',
        'Sunrise', 'Sunset',
        ),
    'ProPlot Diverging': (
        'NegPos', 'Div', 'DryWet', 'Moisture',
        ),
    # Nice diverging maps
    'Miscellaneous Diverging': (
        'ColdHot', 'CoolWarm', 'BR',
        ),
    # cmOcean
    'cmOcean Sequential': (
        'Oxy', 'Thermal', 'Dense', 'Ice', 'Haline',
        'Deep', 'Algae', 'Tempo', 'Speed', 'Turbid', 'Solar', 'Matter',
        'Amp', 'Phase',
        ),
    'cmOcean Diverging': (
        'Balance', 'Delta', 'Curl',
        ),
    # ColorBrewer
    'ColorBrewer2.0 Sequential': (
        'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'PuBu', 'PuBuGn', 'BuGn', 'GnBu', 'YlGnBu', 'YlGn'
        ),
    'ColorBrewer2.0 Diverging': (
        'Spectral', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGY',
        'RdBu', 'RdYlBu', 'RdYlGn',
        ),
    # SciVisColor
    'SciVisColor Blues': (
        'Blue0', 'Blue1', 'Blue2', 'Blue3', 'Blue4', 'Blue5', 'Blue6', 'Blue7', 'Blue8', 'Blue9', 'Blue10', 'Blue11',
        ),
    'SciVisColor Greens': (
        'Green1', 'Green2', 'Green3', 'Green4', 'Green5', 'Green6', 'Green7', 'Green8',
        ),
    'SciVisColor Oranges': (
        'Orange1', 'Orange2', 'Orange3', 'Orange4', 'Orange5', 'Orange6', 'Orange7', 'Orange8',
        ),
    'SciVisColor Browns': (
        'Brown1', 'Brown2', 'Brown3', 'Brown4', 'Brown5', 'Brown6', 'Brown7', 'Brown8', 'Brown9',
        ),
    'SciVisColor Reds/Purples': (
        'RedPurple1', 'RedPurple2', 'RedPurple3', 'RedPurple4', 'RedPurple5', 'RedPurple6', 'RedPurple7', 'RedPurple8',
        ),
    # Builtin maps that will be deleted; categories are taken from comments in
    # matplotlib source code. Some of these are really bad, some are segmented
    # maps when the should be color cycles, and some are just uninspiring.
    'MATLAB': (
        'bone', 'cool', 'copper', 'autumn', 'flag', 'prism',
        'jet', 'hsv', 'hot', 'spring', 'summer', 'winter', 'pink', 'gray',
        ),
    'GNUplot': (
        'gnuplot', 'gnuplot2', 'ocean', 'afmhot', 'rainbow',
        ),
    'GIST': (
        'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar',
        'gist_rainbow', 'gist_stern', 'gist_yarg',
        ),
    'Miscellaneous': (
        'binary', 'bwr', 'brg', # appear to be custom matplotlib, very simple construction
        'cubehelix', 'wistia',  'CMRmap', # individually released
        'seismic', 'terrain', 'nipy_spectral', # origin ambiguous
        ),
    }
CMAPS_DELETE = (
    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
    'spring', 'summer', 'autumn', 'winter', 'cool', 'wistia',
    'afmhot', 'gist_heat', 'copper',
    'seismic', 'bwr', 'brg',
    'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
    'gnuplot', 'gnuplot2', 'cmrmap', 'hsv', 'hot', 'rainbow',
    'gist_rainbow', 'jet', 'nipy_spectral', 'gist_ncar', 'cubehelix',
    )
CMAPS_DIV_SLICES = {
    'piyg': (None, 2, None),
    'prgn': (None, 1, 2, None), # purple red green
    'brbg': (None, 2, 3, None), # brown blue green
    'puor': (None, 2, None),
    'rdgy': (None, 2, None),
    'rdbu': (None, 2, None),
    'rdylbu': (None, 2, 4, None),
    'rdylgn': (None, 2, 4, None),
    'br': (None, 1, None),
    'coldhot': (None, 4, None),
    'negpos': (None, 3, None),
    'drywet': (None, 3, None),
    } # slice args used to split up segments of names
CMAPS_DIV_PAIRS = [
    (name, ''.join(reversed([name[slice(*idxs[i:i+2])] for i in range(len(idxs)-1)])),)
    for name,idxs in CMAPS_DIV_SLICES.items()
    ] # tuple pairs of mirror image cmap names

# Color cycle stuff
CYCLES_PRESET = {
    # Default matplotlib v2
    'default': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    # From stylesheets
    '538': ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'],
    'ggplot': ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
    # The two nice-looking seaborn color cycles
    'ColorBlind': ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9'],
    'ColorBlind10': ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"], # versions with more colors
    # Created with iwanthue and coolers
    'FlatUI': ["#3498db", "#e74c3c", "#95a5a6", "#34495e", "#2ecc71", "#9b59b6"],
    'Warm': [(51,92,103), (158,42,43), (255,243,176), (224,159,62), (84,11,14)],
    'Cool': ["#6C464F", "#9E768F", "#9FA4C4", "#B3CDD1", "#C7F0BD"],
    'Sharp': ["#007EA7", "#D81159", "#B3CDD1", "#FFBC42", "#0496FF"],
    'Hot': ["#0D3B66", "#F95738", "#F4D35E", "#FAF0CA", "#EE964B"],
    'Contrast': ["#2B4162", "#FA9F42", "#E0E0E2", "#A21817", "#0B6E4F"],
    'Floral': ["#23395B", "#D81E5B", "#FFFD98", "#B9E3C6", "#59C9A5"],
    }
CYCLES_DELETE = (
    'tab10', 'tab20', 'tab20b', 'tab20c',
    'paired', 'pastel1', 'pastel2', 'dark2',
    ) # unappealing cycles, and cycles that are just merged monochrome colormaps
CYCLES_RENAME = (
    ('Accent','Set1'),
    ) # rename existing cycles

# Named color filter props
FILTER_SPACE = 'hcl' # dist 'distinct-ness' of colors using this colorspace
FILTER_THRESH = 0.10 # bigger number equals fewer colors
FILTER_TRANSLATIONS = tuple((re.compile(regex), sub) for regex,sub in (
    ('/', ' '), ("'s", ''),
    ('grey', 'gray'),
    ('pinky', 'pink'),
    ('greeny', 'green'),
    ('bluey',  'blue'),
    ('purply', 'purple'),
    ('purpley', 'purple'),
    ('yellowy', 'yellow'),
    ('robin egg', 'robins egg'),
    ('egg blue', 'egg'),
    (r'reddish', 'red'),
    (r'purplish', 'purple'),
    (r'bluish',  'blue'),
    (r'ish\b', ''),
    ('bluegray', 'blue gray'),
    ('grayblue', 'gray blue'),
    ('lightblue', 'light blue')
    )) # prevent registering similar-sounding names
FILTER_ADD = (
    'charcoal', 'sky blue', 'eggshell', 'sea blue', 'coral', 'aqua', 'tomato red', 'brick red', 'crimson',
    'red orange', 'yellow orange', 'yellow green', 'blue green',
    'blue violet', 'red violet',
    ) # common names that should always be included
FILTER_BAD = re.compile('(' + '|'.join((
    'shit', 'poop', 'poo', 'pee', 'piss', 'puke', 'vomit', 'snot', 'booger', 'bile', 'diarrhea',
    )) + ')') # filter these out, let's try to be professional here...

# Named color stuff
OPEN_COLORS = (
    'red', 'pink', 'grape', 'violet',
    'indigo', 'blue', 'cyan', 'teal',
    'green', 'lime', 'yellow', 'orange', 'gray'
    )
BASE_COLORS_FULL = {
    'blue':    (0, 0, 1),
    'green':   (0, 0.5, 0),
    'red':     (1, 0, 0),
    'cyan':    (0, 0.75, 0.75),
    'magenta': (0.75, 0, 0.75),
    'yellow':  (0.75, 0.75, 0),
    'black':   (0, 0, 0),
    'white':   (1, 1, 1),
    }

# Color math and color spaces
CSPACE_ALIASES = {
    'rgb':   'rgb',
    'hsv':   'hsv',
    'hpl':   'hpl',
    'hpluv': 'hpl',
    'hsl':   'hsl',
    'hsluv': 'hsl',
    'hcl':   'hcl',
    }
CSPACE_IDXS = {
    'hue': 0,
    'chroma': 1,
    'saturation': 1,
    'luminance': 2,
    'alpha': 3,
    }

# Check data
DATA_USER_PATHS = {*()}
DATA_ALLOWED_PATHS = {'cmaps', 'cycles', 'fonts'}
def _check_data():
    """Checks the data folder and issues helpful warning message for new users.
    This is called inside every register function."""
    global DATA_USER_PATHS
    data_user = os.path.join(os.path.expanduser('~'), '.proplot')
    data_user_paths = {os.path.basename(path) for path in glob.glob(os.path.join(data_user, '*'))}
    if data_user_paths != DATA_USER_PATHS and data_user_paths > DATA_ALLOWED_PATHS:
        DATA_USER_PATHS = data_user_paths
        warnings.warn(f'Found extra files {", ".join(data_user_paths - DATA_ALLOWED_PATHS)} in the ~/.proplot folder. Files must be placed in the .proplot/cmaps, .proplot/cycles, or .proplot/fonts subfolders.')

#-----------------------------------------------------------------------------#
# Classes
#-----------------------------------------------------------------------------#
# Flexible color names class
# * Matplotlib 'color' arguments are passed to to_rgba, which tries
#   to read directly from cache and if that fails, tries to sanitize input.
#   The sanitization raises error when encounters (colormap, idx) tuple. So
#   we need to override the *cache* instead of color dictionary itself!
# * Builtin to_rgb tries to get cached colors as dict[name, alpha],
#   resulting in key as (colorname, alpha) or ((R,G,B), alpha) tuple. Impossible
#   to differentiate this from (cmapname, index) usage! Must do try except lookup
#   into colormap dictionary every time. Don't want to do this for actual
#   color dict for sake of speed, so we only wrap *cache* lookup. Also we try
#   to avoid cmap lookup attempt whenever possible with if statement.
class ColorCacheDict(dict):
    """Special dictionary that lets user draw single color tuples from
    arbitrary colormaps or color cycles."""
    def __getitem__(self, key):
        """
        Either samples the color from a colormap or color cycle,
        or calls the parent getitem to look up the color name.

        For a **smooth colormap**, usage is e.g.
        ``color=('Blues', 0.8)`` -- the number should be between 0 and 1, and
        indicates where to draw the color from the smooth colormap. For a
        "listed" colormap, i.e. a **color cycle**, usage is e.g.
        ``color=('colorblind', 2)``. The number indicates the index in the
        list of discrete colors.

        These examples work with any matplotlib command that accepts
        a ``color`` keyword arg.
        """
        # Pull out alpha and draw color from cmap
        rgb, alpha = key
        if not isinstance(rgb, str) and np.iterable(rgb) and len(rgb) == 2 and isinstance(rgb[1], Number) and isinstance(rgb[0], str): # i.e. is not None; this is *very common*, so avoids lots of unnecessary lookups!
            try:
                cmap = mcm.cmap_d[rgb[0]]
            except (TypeError, KeyError):
                pass
            else:
                if isinstance(cmap, mcolors.ListedColormap):
                    if not 0 <= rgb[1] < len(cmap.colors):
                        raise ValueError(f'Color cycle sample for {rgb[0]!r} cycle must be between 0 and {len(cmap.colors)-1}, got {rgb[1]}.')
                    rgb = cmap.colors[rgb[1]] # draw color from the list of colors, using index
                else:
                    if not 0 <= rgb[1] <= 1:
                        raise ValueError(f'Colormap sample for {rgb[0]!r} colormap must be between 0 and 1, got {rgb[1]}.')
                    rgb = cmap(rgb[1]) # interpolate color from colormap, using key in range 0-1
                rgba = mcolors.to_rgba(rgb, alpha)
                return rgba
        return super().__getitem__((rgb, alpha))

class _ColorMappingOverride(mcolors._ColorMapping):
    def __init__(self, mapping):
        """Wraps the cache."""
        super().__init__(mapping)
        self.cache = ColorCacheDict({})

# Apply subclass
# Modify colorConverter and use that everywhere in ProPlot, so only have to
# reference private API in these three lines.
if not isinstance(mcolors._colors_full_map, _ColorMappingOverride):
    _map = _ColorMappingOverride(mcolors._colors_full_map)
    mcolors._colors_full_map = _map
    mcolors.colorConverter.cache = _map.cache # re-instantiate
    mcolors.colorConverter.colors = _map # re-instantiate

# Flexible colormap identification
class CmapDict(dict):
    def __init__(self, kwargs):
        """
        Flexible, case-insensitive colormap identification. Replaces the
        `matplotlib.cm.cmap_d` dictionary that stores registered colormaps.

        Behaves like a dictionary, with three new features:

        1. Names are case insensitive: ``'Blues'``, ``'blues'``, and ``'bLuEs'``
           are all valid names for the "Blues" colormap.
        2. "Reversed" colormaps are not stored directly: Requesting e.g.
           ``'Blues_r'`` will just look up ``'Blues'``, then return the result
           of the `~matplotlib.colors.Colormap.reversed` method.
        3. Diverging colormap names can be referenced by their "inverse" name.
           For example, ``'BuRd'`` is equivalent to ``'RdBu_r'``, as are
           ``'BuYlRd'`` and ``'RdYlBu_r'``.
        """
        for key,value in kwargs.items():
            if not isinstance(key, str):
                raise KeyError(f'Invalid key {key}. Must be string.')
            if key[-2:] == '_r': # do not need to store these!
                continue
            self[key] = value

    def __getitem__(self, key):
        """Sanitizes key name then queries the dictionary."""
        key = self._sanitize_key(key, mirror=True)
        return self._getitem(key)

    def __setitem__(self, key, item):
        """Sanitizes key name then assigns to the dictionary."""
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key}. Must be string.')
        key = self._sanitize_key(key, mirror=False)
        return super().__setitem__(key, item)

    def __contains__(self, item):
        """Use sanitized key name for `'in'`."""
        try: # by default __contains__ uses object.__getitem__ and ignores overrides
            self.__getitem__(item)
            return True
        except KeyError:
            return False

    def _getitem(self, key, *args):
        """Gets value but skips key sanitization."""
        reverse = False
        if key[-2:] == '_r':
            key, reverse = key[:-2], True
        try:
            value = super().__getitem__(key) # may raise keyerror
        except KeyError as err:
            if len(args) == 1:
                return args[0]
            else:
                raise err
        if reverse:
            try:
                value = value.reversed()
            except AttributeError:
                raise KeyError(f'Dictionary value in {key} must have reversed() method.')
        return value

    def _sanitize_key(self, key, mirror=True):
        """Sanitizes key name."""
        if not isinstance(key, str):
            raise ValueError(f'Invalid key {key}. Must be string.')
        key = key.lower()
        reverse = False
        if key[-2:] == '_r':
            key = key[:-2]
            reverse = True
        if mirror and not super().__contains__(key): # search for mirrored key
            key_mirror = key
            for mirror in CMAPS_DIV_PAIRS:
                try:
                    idx = mirror.index(key)
                    key_mirror = mirror[1 - idx]
                except (ValueError,KeyError):
                    continue
            if super().__contains__(key_mirror):
                reverse = (not reverse)
                key = key_mirror
        if reverse:
            key = key + '_r'
        return key

    def get(self, key, *args):
        """Retrieves sanitized key name."""
        key = self._sanitize_key(key, mirror=True)
        return super().get(key, *args)

    def pop(self, key, *args):
        """Pops sanitized key name."""
        key = self._sanitize_key(key, mirror=True)
        return super().pop(key, *args)

    def update(self, *args, **kwargs):
        """Replicates dictionary update with sanitized key names."""
        if len(args) == 1:
            kwargs.update(args[0])
        elif len(args) > 1:
            raise TypeError(f'update() expected at most 1 arguments, got {len(args)}.')
        for key,value in kwargs.items():
            self[key] = value

# Apply subclass
if not isinstance(mcm.cmap_d, CmapDict):
    mcm.cmap_d = CmapDict(mcm.cmap_d)

#-----------------------------------------------------------------------------#
# Color manipulation functions
#-----------------------------------------------------------------------------#
def _get_space(space):
    """Verify requested colorspace is valid."""
    space = CSPACE_ALIASES.get(space.lower(), None)
    if space is None:
        raise ValueError(f'Unknown colorspace "{space}".')
    return space

def _get_channel(color, channel, space='hsl'):
    """Gets hue, saturation, or luminance channel value from registered
    string color name. The color name `color` can optionally be a string
    with the format ``'color+x'`` or ``'color-x'``, where `x` specifies
    the offset from the channel value."""
    # Interpret channel
    if callable(color) or isinstance(color, Number):
        return color
    channel = CSPACE_IDXS.get(channel, None)
    if channel is None:
        raise ValueError(f'Unknown channel {channel}.')
    if channel not in (0,1,2):
        raise ValueError('Channel must be in [0,1,2].')
    # Interpret string or RGB tuple
    offset = 0
    if isinstance(color, str):
        match = re.search('([-+][0-9.]+)$', color)
        if match:
            offset = float(match.group(0))
            color = color[:match.start()]
    return offset + to_xyz(to_rgb(color), space)[channel]

def shade(color, scale=0.5):
    """Changes the "shade" of a color by scaling its luminance channel by `scale`."""
    color = to_rgb(color) # ensure is valid color
    color = [*colormath.rgb_to_hsl(*color)]
    color[2] = max(0, min(color[2]*scale, 100)) # multiply luminance by this value
    color = [*colormath.hsl_to_rgb(*color)]
    return tuple(color)

def saturate(color, scale=0.5):
    """Changes the saturation of a color by scaling its saturation channel by `scale`."""
    color = to_rgb(color) # ensure is valid color
    color = [*colormath.rgb_to_hsl(*color)]
    color[1] = max(0, min(color[1]*scale, 100)) # multiply luminance by this value
    color = [*colormath.hsl_to_rgb(*color)]
    return tuple(color)

def to_rgb(color, space='rgb', cycle=None):
    """Generalization of matplotlib's `~matplotlib.colors.to_rgb`. Translates
    colors from *any* colorspace to RGB, converts color strings to RGB
    tuples, and transforms color cycle strings (e.g. ``'C0'``, ``'C1'``, ``'C2'``)
    into their corresponding RGB colors using the input `cycle`, which defaults
    to the current color cycler. Inverse of `to_xyz`."""
    # Convert color cycle strings
    if isinstance(color, str) and re.match('^C[0-9]$', color):
        if isinstance(cycle, str):
            try:
                cycle = mcm.cmap_d[cycle].colors
            except Exception:
                cycles = sorted(name for name,cmap in mcm.cmap_d.items() if isinstance(cmap, mcolors.ListedColormap))
                raise ValueError(f'Invalid cycle name "{cycle}". Options are: {", ".join(cycles)}')
        elif cycle is None:
            cycle = rcParams['axes.prop_cycle'].by_key()
            if 'color' not in cycle:
                cycle = ['k']
            else:
                cycle = cycle['color']
        elif not np.iterable(cycle):
            raise ValueError(f'Invalid cycle "{cycle}".')
        color = cycle[int(color[-1]) % len(cycle)]
    # Translate RGB strings and (cmap,index) tuples
    if isinstance(color, str) or (np.iterable(color) and len(color) == 2):
        try:
            color = mcolors.to_rgb(color) # ensure is valid color
        except Exception:
            raise ValueError(f'Invalid RGB argument "{color}".')
    elif space == 'rgb':
        color = color[:3] # trim alpha
        try:
            if any(c > 1 for c in color):
                color = [c/255 for c in color] # scale to within 0-1
            color = tuple(color)
        except Exception:
            raise ValueError(f'Invalid RGB argument {color}.')
    # Translate from other colorspaces
    elif space == 'hsv':
        color = colormath.hsl_to_rgb(*color)
    elif space == 'hpl':
        color = colormath.hpluv_to_rgb(*color)
    elif space == 'hsl':
        color = colormath.hsluv_to_rgb(*color)
    elif space == 'hcl':
        color = colormath.hcl_to_rgb(*color)
    else:
        raise ValueError('Invalid color "{color}" for colorspace "{space}".')
    return color

def to_xyz(color, space):
    """Translates from the RGB colorspace to colorspace `space`. Inverse
    of `to_rgb`."""
    # Run tuple conversions
    # NOTE: Don't pass color tuple, because we may want to permit out-of-bounds RGB values to invert conversion
    color = to_rgb(color)
    if space == 'hsv':
        color = colormath.rgb_to_hsl(*color) # rgb_to_hsv would also work
    elif space == 'hpl':
        color = colormath.rgb_to_hpluv(*color)
    elif space == 'hsl':
        color = colormath.rgb_to_hsluv(*color)
    elif space == 'hcl':
        color = colormath.rgb_to_hcl(*color)
    elif space == 'rgb':
        pass
    else:
        raise ValueError(f'Invalid colorspace {space}.')
    return color

#-----------------------------------------------------------------------------#
# Helper functions
#-----------------------------------------------------------------------------#
def _clip_colors(colors, clip=True, gray=0.2):
    """
    Clips impossible colors rendered in an HSl-to-RGB colorspace conversion.
    Used by `PerceptuallyUniformColormap`. If `mask` is ``True``, impossible
    colors are masked out

    Parameters
    ----------
    colors : list of length-3 tuples
        The RGB colors.
    clip : bool, optional
        If `clip` is ``True`` (the default), RGB channel values >1 are clipped
        to 1. Otherwise, the color is masked out as gray.
    gray : float, optional
        The identical RGB channel values (gray color) to be used if `mask`
        is ``True``.
    """
    # Clip colors
    colors = np.array(colors)
    over = (colors > 1)
    under = (colors < 0)
    if clip:
        colors[under] = 0
        colors[over]  = 1
    else:
        colors[(under | over)] = gray
    # Message
    # NOTE: Never print warning because happens when using builtin maps
    # message = 'Clipped' if clip else 'Invalid'
    # for i,name in enumerate('rgb'):
    #     if under[:,i].any():
    #         warnings.warn(f'{message} "{name}" channel ( < 0).')
    #     if over[:,i].any():
    #         warnings.warn(f'{message} "{name}" channel ( > 1).')
    return colors

def _slice_cmap(cmap, left=None, right=None, name=None, N=None):
    """Helper function that cleanly divides linear segmented colormaps and
    subsamples listed colormaps. Full documentation is in `Colormap`."""
    # Bail out
    if left is None and right is None:
        return cmap

    # Simple process for listed colormap, just truncate the colors
    name = name or 'no_name'
    if isinstance(cmap, mcolors.ListedColormap):
        try:
            return mcolors.ListedColormap(cmap.colors[left:right])
        except Exception:
            raise ValueError(f'Invalid slice {slice(left,right)} for listed colormap.')

    # Resample the segmentdata arrays
    data = {}
    left = left or 0
    right = right or 1
    dict_ = {key:value for key,value in cmap._segmentdata.items() if 'gamma' not in key}
    gammas = {'saturation':'gamma1', 'luminance':'gamma2'}
    kwargs = {}
    for key,xyy in dict_.items():
        # Get coordinates
        xyy = np.array(xyy)
        x   = xyy[:,0]
        xleft,  = np.where(x > left)
        xright, = np.where(x < right)
        if len(xleft) == 0:
            raise ValueError(f'Invalid x minimum {left}.')
        if len(xright) == 0:
            raise ValueError(f'Invalid x maximum {right}.')
        # Slice
        # l is the first point where x > 0 or x > left, should be at least 1
        # r is the last point where r < 1 or r < right
        l, r = xleft[0], xright[-1]
        ixyy = xyy[l:r+1,:].copy()
        xl = xyy[l-1,1:] + (left - x[l-1])*(xyy[l,1:] - xyy[l-1,1:])/(x[l] - x[l-1])
        ixyy = np.concatenate(([[left, *xl]], ixyy), axis=0)
        xr = xyy[r,1:] + (right - x[r])*(xyy[r+1,1:] - xyy[r,1:])/(x[r+1] - x[r])
        ixyy = np.concatenate((ixyy, [[right, *xr]]), axis=0)
        ixyy[:,0] = (ixyy[:,0] - left)/(right - left)
        data[key] = ixyy
        # Retain the corresponding 'gamma' *segments*
        if key not in gammas:
            continue
        gamma = getattr(cmap, '_' + gammas[key])
        if np.iterable(gamma):
            gamma = gamma[l-1:r+1]
        kwargs[gammas[key]] = gamma

    # Rebuild cmap
    if hasattr(cmap, '_space'):
        kwargs['space'] = cmap._space
    return type(cmap)(name, data, N=cmap.N, **kwargs)

def _shift_cmap(cmap, shift=None, name=None):
    """Shift a cyclic colormap by `shift` degrees out of 360 degrees."""
    # Bail out
    if not shift:
        return cmap

    # Rotate colors for listed colormap
    name = name or 'no_name'
    if isinstance(cmap, mcolors.ListedColormap):
        shift = shift % len(cmap.colors)
        colors = [*cmap.colors] # ensure list
        colors = colors[shift:] + colors[:shift]
        return mcolors.ListedColormap(colors, name=name, N=len(colors))

    # Trickier for smooth colormaps, must shift coordinates
    # TODO: This won't work for lo-res colormaps or percpetually
    # uniform maps with only two coordinate transitions, right?
    data = cmap._segmentdata.copy()
    for key,orig in cmap._segmentdata.items():
        # Drop an end color
        orig = np.array(orig)
        orig = orig[1:,:]
        array = orig.copy()
        array[:,0] -= shift/360
        array[:,0] %= 1
        # Add end color back in
        array = array[array[:,0].argsort(),:]
        array = np.concatenate((array[-1:,:], array), axis=0)
        array[:1,0] = array[1:2,0] - np.diff(array[1:3,0])
        # Normalize x-range
        array[:,0] -= array[:,0].min()
        array[:,0] /= array[:,0].max()
        data[key] = array

    # Generate shifted colormap
    cmap = mcolors.LinearSegmentedColormap(name, data, N=cmap.N)
    cmap._cyclic = True
    return cmap

def _merge_cmaps(*imaps, ratios=1, name=None, N=512, **kwargs):
    """Merges arbitrary colormaps. This is used when you pass multiple `imaps`
    to the `Colormap` function. Full documentation is in `Colormap`."""
    # Bail out
    if len(imaps) == 1:
        return imaps[0]
    types = {type(cmap) for cmap in imaps}
    if len(types) != 1:
        raise ValueError(f'Mixed colormap types {types}. Maps must all be LinearSegmentedColormap or PerceptuallyUniformColormap.')
    type_ = types.pop()

    # For listed colormap, just combine the colors
    name = name or 'no_name'
    if all(isinstance(cmap, mcolors.ListedColormap) for cmap in imaps):
        colors = [color for cmap in imaps for color in cmap.colors]
        return mcolors.ListedColormap(colors, name=name, N=len(colors))

    # Tricker for smooth maps
    kwargs = {}
    segmentdata = {}
    ratios = ratios or 1
    if isinstance(ratios, Number):
        ratios = [1]*len(imaps)
    ratios = np.array(ratios)/np.sum(ratios) # so if 4 cmaps, will be 1/4
    x0 = np.concatenate([[0], np.cumsum(ratios)]) # coordinates for edges
    xw = x0[1:] - x0[:-1] # widths between edges

    # PerceptuallyUniformColormaps checks
    if type_ is PerceptuallyUniformColormap:
        # Colorspace
        spaces = {cmap._space for cmap in imaps}
        if len(spaces) > 1:
            raise ValueError(f'Cannot merge colormaps in the different HSL spaces {repr(spaces)}.')
        kwargs['space'] = spaces.pop()
        # Gamma handling
        gammas = {0:'saturation', 1:'luminance'}
        for i,key in enumerate(('gamma1', 'gamma2')):
            if key not in kwargs:
                kwargs[key] = []
            for cmap in imaps:
                gamma = getattr(cmap, '_' + key)
                if not np.iterable(gamma):
                    gamma = [gamma]*(len(cmap._segmentdata[gammas[i]])-1) # length is *number* of rows in segmentdata
                kwargs[key].extend([*gamma])

    # Combine the segmentdata, and use the y1/y2 slots at merge points so
    # we never interpolate between end colors of different colormaps
    keys = {key for cmap in imaps for key in cmap._segmentdata.keys()}
    for key in keys:
        # Combine xyy data
        # WARNING: If just reference a global 'funcs' list from inside the
        # 'data' function, end up with grayscale colormap because each 'data'
        # function reads 'funcs' as from the final channel in 'keys'. Must
        # embed 'funcs' into each definition using a keyword argument.
        callable_ = [callable(cmap._segmentdata[key]) for cmap in imaps]
        if all(callable_): # expand range from x-to-w to 0-1
            funcs = [cmap._segmentdata[key] for cmap in imaps]
            def data(ix, funcs=funcs):
                ix = np.atleast_1d(ix)
                kx = np.empty(ix.shape)
                for j,jx in enumerate(ix.flat):
                    idx = max(np.searchsorted(x0, jx)-1, 0)
                    kx.flat[j] = funcs[idx]((jx - x0[idx])/xw[idx])
                return kx
        elif not any(callable_):
            datas = []
            for x,w,cmap in zip(x0[:-1], xw, imaps):
                data = np.array(cmap._segmentdata[key])
                data[:,0] = x + w*data[:,0]
                datas.append(data)
            for i in range(len(datas)-1):
                datas[i][-1,2] = datas[i+1][0,2] # jump to next colormap, never interpolate between colors from different maps
                datas[i+1] = datas[i+1][1:,:] # shave off initial color on next colormap
            data = np.concatenate(datas, axis=0)
            data[:,0] = data[:,0]/data[:,0].max(axis=0) # scale to make maximum exactly 1 (avoid floating point errors)
        else:
            raise ValueError('Mixed callable and non-callable colormap values.')
        segmentdata[key] = data
    return type_(name, segmentdata, N=N, **kwargs)

def _make_segmentdata_array(values, ratios=None, reverse=False, **kwargs):
    """Constructs a list of linear segments for an individual channel.
    This was made so that user can input e.g. a callable function for
    one channel, but request linear interpolation for another one."""
    # Allow callables
    if callable(values):
        if reverse:
            values = lambda x: values(1-x)
        return values # just return the callable
    values = np.atleast_1d(values)
    if len(values) == 1:
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

def make_mapping_array(N, data, gamma=1.0, reverse=False):
    r"""
    Mostly a copy of `~matplotlib.colors.makeMappingArray`, but allows
    *circular* hue gradations along 0-360, disables clipping of
    out-of-bounds channel values, and with fancier "gamma" scaling.

    Parameters
    ----------
    N : int
        Number of points in the generated lookup table.
    data : 2D array-like
        List of :math:`(x, y_0, y_1)` tuples specifying the channel jump (from
        :math:`y_0` to :math:`y_1`) and the :math:`x` coordinate of that
        transition (ranges between 0 and 1).
        See `~matplotlib.colors.LinearSegmentedColormap` for details.
    gamma : float or list of float, optional
        To obtain channel values between coordinates :math:`x_i` and
        :math:`x_{i+1}` in rows :math:`i` and :math:`i+1` of `data`,
        we use the formula:

        .. math::

            y = y_{1,i} + w_i^{\gamma_i}*(y_{0,i+1} - y_{1,i})

        where :math:`\gamma_i` corresponds to `gamma` and the weight
        :math:`w_i` ranges from 0 to 1 between rows ``i`` and ``i+1``.
        If `gamma` is float, it applies to every transition. Otherwise,
        its length must equal ``data.shape[0]-1``.
    reverse : bool, optional
        If ``True``, :math:`w_i^{\gamma_i}` is replaced with
        :math:`1 - (1 - w_i)^{\gamma_i}` -- that is, when `gamma` is greater
        than 1, this weights colors toward *higher* channel values instead
        of lower channel values.

        This is implemented in case we want to apply *equal* "gamma scaling"
        to different HSL channels in different directions. Usually, this
        is done to weight low data values with higher luminance *and* lower
        saturation, thereby emphasizing "extreme" data values with stronger
        colors.
    """
    # Allow for *callable* instead of linearly interpolating between segments
    gammas = np.atleast_1d(gamma)
    if (gammas < 0.01).any() or (gammas > 10).any():
        raise ValueError('Gamma can only be in range [0.01,10].')
    if callable(data):
        if len(gammas) > 1:
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
    if len(gammas) != 1 and len(gammas) != shape[0]-1:
        raise ValueError(f'Need {shape[0]-1} gammas for {shape[0]}-level mapping array, but got {len(gamma)}.')
    if len(gammas) == 1:
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
    # The ui are starting-points, the ci are counts from that point
    # over which segment applies (i.e. where to apply the gamma), the relevant
    # 'segment' is to the *left* of index returned by searchsorted
    _, uind, cind = np.unique(ind, return_index=True, return_counts=True)
    for ui,ci in zip(uind,cind): # length should be N-1
        gamma = gammas[ind[ui]-1] # the relevant segment is to *left* of this number
        if gamma == 1:
            continue
        ireverse = False
        if ci > 1: # i.e. more than 1 color in this 'segment'
            ireverse = ((y0[ind[ui]] - y1[ind[ui]-1]) < 0) # by default want to weight toward a *lower* channel value
        if reverse:
            ireverse = (not ireverse)
        if ireverse:
            distance[ui:ui + ci] = 1 - (1 - distance[ui:ui + ci])**gamma
        else:
            distance[ui:ui + ci] **= gamma

    # Perform successive linear interpolations all rolled up into one equation
    lut = np.zeros((N,), float)
    lut[1:-1] = distance*(y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0]  = y1[0]
    lut[-1] = y0[-1]
    return lut

#-----------------------------------------------------------------------------#
# Generalized colormap/cycle constructors
#-----------------------------------------------------------------------------#
def colors(*args, **kwargs):
    """Identical to `Cycle`, but returns a list of colors instead of
    a `~cycler.Cycler` object."""
    cycle = Cycle(*args, **kwargs)
    return [dict_['color'] for dict_ in cycle]

def Colormap(*args, name=None, cyclic=None, listmode='perceptual',
        fade=None, cycle=None,
        shift=None, cut=None, left=None, right=None, reverse=False,
        ratios=1, gamma=None, gamma1=None, gamma2=None,
        save=False, N=None,
        **kwargs):
    """
    Function for generating and merging colormaps in a variety of ways;
    used to interpret the `cmap` and `cmap_kw` arguments when passed to
    any plotting method wrapped by `~proplot.wrappers.cmap_wrapper`.

    Parameters
    ----------
    *args : colormap-spec
        Positional arguments that individually generate colormaps. If more than
        one argument is passed, the resulting colormaps are merged. Arguments
        are interpreted as follows.

        * If `~matplotlib.colors.LinearSegmentedColormap`, nothing more is done.
        * If string colormap name, that `~matplotlib.colors.LinearSegmentedColormap`
          is looked up and used.
        * If RGB tuple or color string, a `PerceptuallyUniformColormap` is
          generated with `monochrome_cmap`. If the string ends in ``'_r'``,
          the monochrome map will be *reversed*, i.e. will go from dark to light
          instead of light to dark.
        * If list of RGB tuples or color strings, a `PerceptuallyUniformColormap`
          is generated by passing the list to the `~PerceptuallyUniformColormap.from_list`
          method.
        * If dictionary, a `PerceptuallyUniformColormap` is generated by
          passing the dictionary as keyword args to the
          `~PerceptuallyUniformColormap.from_hsl` method.

    name : str, optional
        Name of the resulting colormap. Default name is ``'no_name'``.
        The resulting colormap can then be invoked by passing ``cmap='name'``
        to plotting functions like `~matplotlib.axes.Axes.contourf`.
    cyclic : bool, optional
        Whether the colormap is cyclic. Will cause `~proplot.wrappers.cmap_wrapper`
        to pass this flag to `BinNorm`. This will prevent having the same color
        on either end of the colormap.
    listmode : {'perceptual', 'linear', 'listed'}, optional
        Controls how colormaps are generated when you input list(s) of colors.
        If ``'perceptual'``, a `PerceptuallyUniformColormap` is generated with
        `PerceptuallyUniformColormap.from_list`. If ``'linear'``,
        a `~matplotlib.colors.LinearSegmentedColormap` is generated with
        `~matplotlib.colors.LinearSegmentedColormap.from_list`. If
        ``'listed'``, the `~matplotlib.colors.ListedColormap` is generated.
        Default is ``'perceptual'`` when calling `Colormap` directly, and
        ``'listed'`` when `Colormap` is called by `Cycle`.
    perceptual : bool, optional
    fade : float, optional
        The maximum luminosity used when generating `monochrome_cmap` colormaps.
        Default is ``100`` when calling `Colormap` directly, and ``90`` when
        `Colormap` is called by `Cycle` (this prevents having pure white in
        the color cycle).

        For example, ``plot.Colormap('blue', fade=80)`` generates a blue
        colormap that fades to a pale blue with 80% luminance.
    cycle : str or list of color-spec, optional
        The registered cycle name or a list of colors used to interpret cycle
        color strings like ``'C0'`` and ``'C2'`` when generating
        `monochrome_cmap` colormaps. Default is colors from the currently
        active property cycler.

        For example, ``plot.Colormap('C0', 'C1', 'C2', cycle='538')``
        generates a colormap using colors from the ``'538'`` color cycle.
    shift : float, optional
        For `~matplotlib.colors.LinearSegmentedColormap` maps, this
        rotates the colors by `shift` degrees out of 360 degrees. This is
        mainly useful for "cyclic" colormaps. For example, ``shift=180``
        moves the edge colors to the center of the colormap.

        For `~matplotlib.colors.ListedColormap` maps, this rotates
        the color list by `shift` places. For example, ``shift=2`` moves the
        start of the color cycle two places to the right.
    left, right : float or list of float, optional
        For `~matplotlib.colors.LinearSegmentedColormap` maps, this
        deletes colors on the left and right sides of the
        colormap(s). For example, ``left=0.1`` deletes the leftmost 10% of
        the colormap, while ``right=0.9`` deletes the rightmost 10%.

        For `~matplotlib.colors.ListedColormap` maps, this slices
        the color list using ``cmap.colors = cmap.colors[left:right]``.
        For example, ``left=1`` with ``right=None`` deletes the first color.

        If float, these apply to the final, *merged* colormap. If list of float,
        these apply to *each* individual colormap before the colormaps are
        merged. There is no difference if ``len(args) == 1``.
    cut : float, optional
        For `~matplotlib.colors.LinearSegmentedColormap` maps, this
        cuts out colors in the *center* of the colormap. This is useful
        if you want to have a sharper cutoff between "negative" and "positive"
        values in a diverging colormap. For example, ``cut=0.1`` cuts out
        the middle 10% of the colormap.
    reverse : bool or list of bool, optional
        Optionally reverses the colormap(s).

        If bool, this applies to the final, *merged* colormap. If list of bool,
        these apply to *each* individual colormap before the colormaps are
        merged. There is no difference if ``len(args) == 1``.
    ratios : list of float, optional
        Indicates the ratios used to *merge* the colormaps. Length must
        equal ``len(args)``. For example, if `args` contains
        ``['blues', 'reds']`` and `ratios` is ``[2, 1]``, this generates a
        colormap with two-thirds blue colors on the left and one-third red
        colors on the right.
    gamma1, gamma2, gamma : float, optional
        Gamma-scaling for the saturation, luminance, and both channels
        for perceptualy uniform colormaps. See the
        `PerceptuallyUniformColormap` documentation.
    save : bool, optional
        Whether to save the colormap in the folder ``~/.proplot``. The
        folder is created if it does not already exist.

        If the colormap is a `~matplotlib.colors.ListedColormap` (i.e. a
        "color cycle"), the list of hex strings are written to
        ``cycles/name.hex``.

        If the colormap is a `~matplotlib.colors.LinearSegmentedColormap`,
        the segment data dictionary is written to ``cmaps/name.json``.
    N : int, optional
        Number of colors to generate in the hidden lookup table ``_lut``.
        By default, a relatively high resolution of 256 is chosen (see notes).

    Returns
    -------
    `~matplotlib.colors.Colormap`
        A `~matplotlib.colors.LinearSegmentedColormap` or
        `~matplotlib.colors.ListedColormap` instance.
    """
    # Initial stuff
    if not args:
        raise ValueError(f'Colormap requires at least one positional argument.')
    if listmode not in ('listed', 'linear', 'perceptual'):
        raise ValueError(f'Invalid listmode={listmode!r}. Options are "listed", "linear", and "perceptual".')
    N_ = N or rcParams['image.lut']
    name = name or 'no_name' # must have name, mcolors utilities expect this
    imaps = []
    for i,cmap in enumerate(args):
        # Retrieve Colormap instance. Makes sure lookup table is reset.
        ireverse = False
        if np.iterable(reverse) and reverse[i]:
            ireverse = True
        if isinstance(cmap,str) and cmap in mcm.cmap_d:
            cmap = mcm.cmap_d[cmap]
            if kwargs:
                warnings.warn(f'Ignoring extra kwargs {kwargs}.')
        if isinstance(cmap, mcolors.ListedColormap):
            if kwargs:
                warnings.warn(f'Ignoring extra kwargs {kwargs}.')
        elif isinstance(cmap, mcolors.LinearSegmentedColormap):
            if kwargs:
                warnings.warn(f'Ignoring extra kwargs {kwargs}.')
            # Resample, allow overriding the gamma and copy over add-on attribute
            cyclic = getattr(cmap, '_cyclic', False)
            cmap = cmap._resample(N_) # makes a copy!
            cmap._cyclic = cyclic
            if isinstance(cmap, PerceptuallyUniformColormap):
                cmap._gamma1 = _notNone(gamma1, gamma, cmap._gamma1)
                cmap._gamma2 = _notNone(gamma2, gamma, cmap._gamma2)
            elif gamma:
                cmap._gamma = _notNone(gamma, cmap._gamma)
        # Build colormap on-the-fly
        elif isinstance(cmap, dict):
            # Dictionary of hue/sat/luminance values or 2-tuples representing linear transition
            try:
                cmap = PerceptuallyUniformColormap.from_hsl(name, N=N_, **cmap, **kwargs)
            except Exception:
                raise ValueError(f'Invalid input "{cmap}" for from_hsl method.')
        elif not isinstance(cmap, str) and np.iterable(cmap) and all(np.iterable(color) for color in cmap):
            # List of color tuples or color strings, i.e. iterable of iterables
            cmap = [to_rgb(color, cycle=cycle) for color in cmap] # transform C0, C1, etc. to actual names
            if listmode == 'listed':
                cmap = mcolors.ListedColormap(cmap, name=name, **kwargs)
            elif listmode == 'linear':
                cmap = mcolors.LinearSegmentedColormap.from_list(name, cmap, **kwargs)
            else:
                cmap = PerceptuallyUniformColormap.from_list(name, cmap, **kwargs)
        else:
            # Monochrome colormap from input color
            if isinstance(cmap, str) and cmap[-2:] == '_r':
                cmap = cmap[:-2]
                ireverse = (not ireverse)
            try:
                color = to_rgb(cmap, cycle=cycle)
            except Exception:
                if not isinstance(cmap, str):
                    raise ValueError(f'Invalid cmap, cycle, or color "{cmap}".')
                else:
                    raise ValueError(f'Invalid cmap, cycle, or color "{cmap}".\n'
                        f'VALID CMAP AND CYCLE NAMES: {", ".join(sorted(mcm.cmap_d))}.\n'
                        f'VALID COLOR NAMES: {", ".join(sorted(mcolors.colorConverter.colors.keys()))}.')
            fade = _notNone(fade, 100)
            cmap = monochrome_cmap(color, fade, name=name, N=N_, **kwargs)

        # Optionally transform colormap by clipping colors or reversing
        if ireverse:
            cmap = cmap.reversed()
        cmap = _slice_cmap(cmap, None if not np.iterable(left) else left[i],
                                None if not np.iterable(right) else right[i], N=N)
        imaps += [cmap]

    # Now merge the result of this arbitrary user input
    # Since we are merging cmaps, potentially *many* color transitions; use big number by default
    N_ = N_*len(imaps)
    if len(imaps) > 1:
        cmap = _merge_cmaps(*imaps, name=name, ratios=ratios, N=N_)

    # Cut out either edge
    left = None if np.iterable(left) else left
    right = None if np.iterable(right) else right
    if not cut: # non-zero and not None
        cmap = _slice_cmap(cmap, left, right, name=name, N=N)
    # Cut out middle colors of a diverging map
    else:
        cright, cleft = 0.5 - cut/2, 0.5 + cut/2
        lcmap = _slice_cmap(cmap, left, cright)
        rcmap = _slice_cmap(cmap, cleft, right)
        cmap = _merge_cmaps(lcmap, rcmap, name=name, N=N_)
    # Cyclic colormap settings
    if shift: # i.e. is non-zero
        cmap = _shift_cmap(cmap, shift, name=name)
    if cyclic is not None:
        cmap._cyclic = cyclic
    elif not hasattr(cmap, '_cyclic'):
        cmap._cyclic = False
    # Optionally reverse
    if not np.iterable(reverse) and reverse:
        cmap = cmap.reversed()
    # Initialize (the _resample methods generate new colormaps,
    # so current one is uninitializied)
    if not cmap._isinit:
        cmap._init()

    # Register the colormap
    mcm.cmap_d[name] = cmap
    # Optionally save colormap to disk
    if save:
        # Save listed colormap i.e. color cycle
        if isinstance(cmap, mcolors.ListedColormap):
            basename = f'{name}.hex'
            filename = os.path.join(DATA_USER_CYCLES, basename)
            with open(filename, 'w') as f:
                f.write(','.join(mcolors.to_hex(color) for color in cmap.colors))
        # Save segment data directly
        else:
            basename = f'{name}.json'
            filename = os.path.join(DATA_USER_CMAPS, basename)
            data = {}
            for key,value in cmap._segmentdata.items():
                data[key] = np.array(value).astype(float).tolist() # from np.float to builtin float, and to list of lists
            for key in ('space','gamma1','gamma2'):
                if hasattr(cmap, '_' + key):
                    data[key] = getattr(cmap, '_' + key)
            with open(filename, 'w') as file:
                json.dump(data, file, indent=4)
        print(f'Saved colormap to "{basename}".')
    return cmap

def Cycle(*args, samples=None, name=None, save=False,
    marker=None, alpha=None, dashes=None, linestyle=None, linewidth=None,
    markersize=None, markeredgewidth=None, markeredgecolor=None, markerfacecolor=None,
    **kwargs):
    """
    Function for generating and merging `~cycler.Cycler` instances in a variety of ways;
    used to interpret the `cycle` and `cycle_kw` arguments when passed to
    any plotting method wrapped by `~proplot.wrappers.cycle_wrapper`.

    If you just want a list of colors instead of a `~cycler.Cycler` instance,
    use the `colors` function. If you want a `~cycler.Cycler` instance that
    imposes black as the default color and cycles through properties like
    ``linestyle`` instead, call this function without any positional arguments.

    Parameters
    ----------
    *args : colormap-spec or cycle-spec, optional
        Positional arguments that control the color values in the `~cycler.Cycler`
        object. If more than one argument is passed, the resulting cycles are
        merged. Arguments are interpreted as follows.

        * If `~cycler.Cycler`, nothing more is done.
        * If list of RGB tuples or color strings, these colors are used.
        * If `~matplotlib.colors.ListedColormap`, colors from the ``colors``
          attribute are used.
        * If string color cycle name, that `~matplotlib.colors.ListedColormap`
          is looked up and its ``colors`` attribute is used.
        * Otherwise, the argument is passed to `Colormap`, and colors
          from the resulting `~matplotlib.colors.LinearSegmentedColormap`
          are used. See the `samples` argument.

        If the last positional argument is numeric, it is used for the `samples`
        keyword argument.
    samples : float or list of float, optional
        For `~matplotlib.colors.ListedColormap` maps, this is the number of
        colors to select. For example, ``Cycle('538', 4)`` returns the first 4
        colors of the ``'538'`` color cycle.

        For `~matplotlib.colors.LinearSegmentedColormap` maps, this is either
        a list of sample coordinates used to draw colors from the map, or an
        integer number of colors to draw. If the latter, the sample coordinates
        are ``np.linspace(0, 1, samples)``. For example, ``Cycle('Reds', 5)``
        divides the ``'Reds'`` colormap into five evenly spaced colors.
    name : str, optional
        Name of the resulting `~matplotlib.colors.ListedColormap` used to
        register the color cycle. Default name is ``'no_name'``.
    save : bool, optional
        Whether to save the color cycle in the folder ``~/.proplot``. The
        folder is created if it does not already exist.

        The cycle is saved as a list of hex strings to the file ``name.hex``.
    marker, alpha, dashes, linestyle, linewidth, markersize, markeredgewidth, markeredgecolor, markerfacecolor : list of specs, optional
        Lists of `~matplotlib.lines.Line2D` properties that can be
        added to the `~cycler.Cycler` instance. If the lists have unequal length,
        they will be filled to match the length of the longest list.
        See `~matplotlib.axes.Axes.set_prop_cycle` for more info on property cyclers.
        Also see the `line style reference <https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html>`__,
        `marker reference <https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/marker_reference.html>`__,
        and the `custom dashes reference <https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/line_demo_dash_control.html>`__.
    **kwargs
        Passed to `Colormap`.

    Returns
    -------
    `~cycler.Cycler`
        A cycler instance that can be passed to `~matplotlib.axes.Axes.set_prop_cycle`.
    """
    # Add properties
    props = {}
    nprops = 0
    for key,value in (
        ('marker',marker),
        ('alpha',alpha),
        ('dashes',dashes),
        ('linestyle',linestyle),
        ('linewidth',linewidth),
        ('markersize',markersize),
        ('markeredgewidth',markeredgewidth),
        ('markeredgecolor',markeredgecolor),
        ('markerfacecolor',markerfacecolor),
        ):
        if value is not None:
            if isinstance(value, str) or not np.iterable(value):
                raise ValueError(f'Invalid {key} property {value}. Must be list or tuple of properties.')
            nprops = max(nprops, len(value))
            props[key] = [*value] # ensure mutable list
    # If args is non-empty, means we want color cycle; otherwise is always black
    if not args:
        props['color'] = ['k'] # ensures property cycler is non empty
        if kwargs:
            warnings.warn(f'Ignoring keyword args {kwargs}.')
    elif all(isinstance(arg, cycler.Cycler) for arg in args):
        # Merge cycler objects
        if kwargs:
            warnings.warn(f'Ignoring keyword args {kwargs}.')
        if len(args) == 1:
            return args[0]
        else:
            props = {}
            for arg in args:
                for key,value in arg.by_key():
                    if key not in props:
                        props[key] = []
                    props[key].extend([*value])
            return cycler.cycler(**props)
    else:
        # Construct and register ListedColormap
        if args and isinstance(args[-1], Number):
            args, samples = args[:-1], args[-1] # means we want to sample existing colormaps or cycles
        kwargs.setdefault('fade', 90)
        kwargs.setdefault('listmode', 'listed')
        cmap = Colormap(*args, **kwargs) # the cmap object itself
        if isinstance(cmap, mcolors.ListedColormap):
            N = samples
            colors = cmap.colors[:N] # if samples is None, does nothing
        else:
            samples = _notNone(samples, 10)
            if isinstance(samples, Integral):
                samples = np.linspace(0, 1, samples) # from edge to edge
            elif np.iterable(samples) and all(isinstance(item,Number) for item in samples):
                samples = np.array(samples)
            else:
                raise ValueError(f'Invalid samples "{samples}".')
            N = len(samples)
            colors = cmap(samples)
        # Register the colormap
        name = name or 'no_name'
        cmap = mcolors.ListedColormap(colors, name=name, N=N)
        cmap.colors = [tuple(color) if not isinstance(color,str) else color for color in cmap.colors] # sanitize
        mcm.cmap_d[name] = cmap
        # Save the cycle
        if save:
            basename = f'{name}.hex'
            filename = os.path.join(DATA_USER_CYCLES, basename)
            with open(filename, 'w') as f:
                f.write(','.join(mcolors.to_hex(color) for color in cmap.colors))
            print(f'Saved color cycle to "{basename}".')
        # Add to property dict
        nprops = max(nprops, len(colors))
        props['color'] = cmap.colors # save the tupled version!
    # Build cycler, make sure lengths are the same
    for key,value in props.items():
        if len(value) < nprops:
            value[:] = [value[i%len(value)] for i in range(nprops)] # make loop double back
    return cycler.cycler(**props)

class PerceptuallyUniformColormap(mcolors.LinearSegmentedColormap):
    """Similar to `~matplotlib.colors.LinearSegmentedColormap`, but instead
    of varying the RGB channels, we vary hue, saturation, and luminance in
    either the HCL colorspace or the HSLuv or HPLuv scalings of HCL."""
    def __init__(self,
        name, segmentdata, space=None, clip=True,
        gamma=None, gamma1=None, gamma2=None, **kwargs
        ):
        """
        Parameters
        ----------
        name : str
            The colormap name.
        segmentdata : dict-like
            Mapping containing the keys ``'hue'``, ``'saturation'``, and
            ``'luminance'``. The key values should be lists containing any of
            the following channel specifiers:

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex tags, in which case the channel
               value for that color is looked up.

            See `~matplotlib.colors.LinearSegmentedColormap` for details.
        space : {'hcl', 'hsl', 'hpl'}, optional
            The hue, saturation, luminance-style colorspace to use for
            interpreting the channels. See `this page
            <http://www.hsluv.org/comparison/>`_ for a description.
        clip : bool, optional
            When we interpolate across HCL space, we can end up with
            "impossible" RGB colors (i.e. RGB channel values >1).

            If `clip` is ``True`` (the default), channel values >1 are clipped
            to 1. Otherwise, the color is masked out as gray.
        gamma1 : float, optional
            If >1, makes low saturation colors more prominent. If <1,
            makes high saturation colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.
        gamma2 : float, optional
            If >1, makes high luminance colors more prominent. If <1,
            makes low luminance colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.
        gamma : float, optional
            Use this to identically set `gamma1` and `gamma2` at once.

        Example
        -------
        The following generates a `PerceptuallyUniformColormap` from a
        `segmentdata` dictionary that uses color names for the hue data,
        instead of channel values between ``0`` and ``360``.

        >>> import proplot as plot
        >>> data = {
        ...     'hue': [[0, 'red', 'red'], [1, 'blue', 'blue']],
        ...     'saturation': [[0, 100, 100], [1, 100, 100]],
        ...     'luminance': [[0, 100, 100], [1, 20, 20]],
        ...     }
        >>> cmap = plot.PerceptuallyUniformColormap(data)

        """
        # Checks
        space = _get_space(_notNone(space, 'hsl'))
        self._space = space
        self._clip  = clip
        keys = {*segmentdata.keys()}
        target = {'hue', 'saturation', 'luminance', 'alpha'}
        if not keys <= target:
            raise ValueError(f'Invalid segmentdata dictionary with keys {keys}.')
        # Gamma scaling
        self._gamma1 = _notNone(gamma1, gamma, 1.0)
        self._gamma2 = _notNone(gamma2, gamma, 1.0)
        # Sanitize segmentdata, convert color strings to their channel values
        for key,array in segmentdata.items():
            if callable(array): # permit callable
                continue
            for i,xyy in enumerate(array):
                xyy = list(xyy) # make copy!
                for j,y in enumerate(xyy[1:]): # modify the y values
                    xyy[j+1] = _get_channel(y, key, space)
                segmentdata[key][i] = xyy
        # Initialize
        super().__init__(name, segmentdata, gamma=1.0, **kwargs)

    def _init(self):
        """As with `~matplotlib.colors.LinearSegmentedColormap`, but converts
        each value in the lookup table from 'input' to RGB."""
        # First generate the lookup table
        channels = ('hue','saturation','luminance')
        reverse = (False, False, True) # gamma weights *low chroma* and *high luminance*
        gammas = (1.0, self._gamma1, self._gamma2)
        self._lut_hsl = np.ones((self.N+3, 4), float) # fill
        for i,(channel,gamma,reverse) in enumerate(zip(channels,gammas,reverse)):
            self._lut_hsl[:-3,i] = make_mapping_array(self.N, self._segmentdata[channel], gamma, reverse)
        if 'alpha' in self._segmentdata:
            self._lut_hsl[:-3,3] = make_mapping_array(self.N, self._segmentdata['alpha'])
        self._lut_hsl[:-3,0] %= 360
        # self._lut_hsl[:-3,0] %= 359 # wrong
        # Make hues circular, set extremes (i.e. copy HSL values)
        self._lut = self._lut_hsl.copy() # preserve this, might want to check it out
        self._set_extremes() # generally just used end values in segmentdata
        self._isinit = True
        # Now convert values to RGBA, and clip colors
        for i in range(self.N+3):
            self._lut[i,:3] = to_rgb(self._lut[i,:3], self._space)
        self._lut[:,:3] = _clip_colors(self._lut[:,:3], self._clip)

    def _resample(self, N):
        """Returns a new colormap with *N* entries."""
        return PerceptuallyUniformColormap(
            self.name, self._segmentdata, space=self._space,
            clip=self._clip, gamma1=self._gamma1, gamma2=self._gamma2, N=N
            )

    def reversed(self, name=None):
        """Returns reversed colormap."""
        if name is None:
            name = self.name + '_r'
        def factory(dat):
            def func_r(x):
                return dat(1.0 - x)
            return func_r
        # Reverse segment data
        data_r = {}
        for key,xyy in self._segmentdata.items():
            if callable(xyy):
                data_r[key] = factory(xyy)
            else:
                data_r[key] = [[1.0 - x, y1, y0] for x, y0, y1 in reversed(xyy)]
        # Reverse gamma scalings, which can be specified on per-segment basis
        kw = {'space':self._space}
        for gamma in ('gamma1','gamma2'):
            kw[gamma] = np.atleast_1d(getattr(self, '_' + gamma))[::-1]
        return PerceptuallyUniformColormap(name, data_r, **kw)

    @staticmethod
    def from_hsl(name, hue=0, saturation=100, luminance=(100, 20), alpha=None, ratios=None, reverse=False, **kwargs):
        """
        Makes a `~PerceptuallyUniformColormap` by specifying the hue, saturation,
        and luminance transitions individually.

        Parameters
        ----------
        hue : float, str, or list thereof, optional
            Hue channel value or list of values. Values can be
            any of the following:

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex strings, in which case the channel
               value for that color is looked up.

            If scalar, the hue does not change across the colormap.
        saturation : float, str, or list thereof, optional
            As with `hue`, but for the saturation channel.
        luminance : float, str, or list thereof, optional
            As with `hue`, but for the luminance channel.
        alpha : float, str, or list thereof, optional
            As with `hue`, but for the alpha channel (the opacity).
        ratios : list of float, optional
            Relative extent of the transitions indicated by the channel
            value lists.

            For example, ``luminance=[100,50,0]`` with ``ratios=[2,1]``
            places the *x*-coordinate where the luminance is 50 at 0.66 --
            the white to gray transition is "slower" than the gray to black
            transition.
        reverse : bool, optional
            Whether to reverse the final colormap.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.
        """
        # Build dictionary, easy peasy
        cdict = {}
        alpha = _notNone(alpha, 1.0)
        for key,channel in zip(('hue','saturation','luminance','alpha'), (hue,saturation,luminance,alpha)):
            cdict[key] = _make_segmentdata_array(channel, ratios, reverse, **kwargs)
        return PerceptuallyUniformColormap(name, cdict, **kwargs)

    @staticmethod
    def from_list(name, colors, ratios=None, reverse=False, **kwargs):
        """
        Makes a `PerceptuallyUniformColormap` from a list of RGB colors.

        Parameters
        ----------
        name : str
            The colormap name.
        colors : list of color-spec
            The list of RGB colors, HEX strings, or registered color names.
        ratios : list of float, optional
            Length ``len(colors)-1`` list of scales for *x*-coordinate
            transitions between colors. Bigger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.
        reverse : bool, optional
            Whether to reverse the result.
        """
        # Translate colors
        # TODO: Allow alpha
        cdict = {}
        space = kwargs.get('space', 'hsl') # use the builtin default
        colors = [to_xyz(color, space) for color in colors]
        channels = [*zip(*colors)]
        if len(channels) not in (3,4):
            raise ValueError(f'Bad color list: {colors}')
        keys = ['hue', 'saturation', 'luminance']
        if len(channels) == 4:
            keys += ['alpha']
        else:
            cdict['alpha'] = lambda x: 1.0 # dummy function that always returns 1.0
        # Build data arrays
        for key,channel in zip(keys,channels):
            cdict[key] = _make_segmentdata_array(channel, ratios, reverse, **kwargs)
        return PerceptuallyUniformColormap(name, cdict, **kwargs)

def monochrome_cmap(color, fade, reverse=False, space='hsl', name='no_name', **kwargs):
    """
    Makes a monochromatic "sequential" colormap that blends from near-white
    to the input color.

    Parameters
    ----------
    color : str or (R,G,B) tuple
        Color RGB tuple, hex string, or named color string.
    fade : float or str or (R,G,B) tuple
        The luminance channel strength, or color from which to take the final
        luminance and saturation. If the former is provided, the saturation
        will be held constant throughout the colormap.
    reverse : bool, optional
        Whether to reverse the colormap.
    space : {'hcl', 'hsl', 'hpl'}, optional
        Colorspace in which the luminance is varied.
    name : str, optional
        Colormap name. Default is ``'no_name'``.

    Other parameters
    ----------------
    **kwargs
        Passed to `PerceptuallyUniformColormap.from_hsl` static method.
    """
    # Get colorspace and channel values
    h, s, l = to_xyz(to_rgb(color), space)
    if isinstance(fade, Number):
        fs, fl = s, fade # fade to *same* saturation by default
    else:
        _, fs, fl = to_xyz(to_rgb(fade), space) # fade to this saturation and this luminance
    # Build colormap
    if reverse:
        s, l = (s,fs), (l,fl) # from color to faded
    else:
        s, l = (fs,s), (fl,l) # from faded to color
    return PerceptuallyUniformColormap.from_hsl(name, h, s, l, space=space, **kwargs)

#-----------------------------------------------------------------------------#
# Return arbitrary normalizer
#-----------------------------------------------------------------------------#
def Norm(norm, levels=None, **kwargs):
    """
    Returns an arbitrary `~matplotlib.colors.Normalize` instance, used to
    interpret the `norm` and `norm_kw` arguments when passed to any plotting
    method wrapped by `~proplot.wrappers.cmap_wrapper`.

    Parameters
    ----------
    norm : str or `~matplotlib.colors.Normalize`
        Key name for the normalizer. The recognized normalizer key names
        are as follows.

        ===============================  ===============================
        Key(s)                           Class
        ===============================  ===============================
        ``'midpoint'``, ``'zero'``       `MidpointNorm`
        ``'segments'``, ``'segmented'``  `LinearSegmentedNorm`
        ``'none'``, ``'null'``           `~matplotlib.colors.NoNorm`
        ``'linear'``                     `~matplotlib.colors.Normalize`
        ``'log'``                        `~matplotlib.colors.LogNorm`
        ``'power'``                      `~matplotlib.colors.PowerNorm`
        ``'symlog'``                     `~matplotlib.colors.SymLogNorm`
        ===============================  ===============================

    levels : array-like, optional
        Level *edges*, passed to `LinearSegmentedNorm` or used to determine
        the `vmin` and `vmax` arguments for `MidpointNorm`.
    **kwargs
        Passed to the `~matplotlib.colors.Normalize` initializer.
        See `this tutorial <https://matplotlib.org/tutorials/colors/colormapnorms.html>`_
        for more info.

    Returns
    -------
    `~matplotlib.colors.Normalize`
        A `~matplotlib.colors.Normalize` instance.
    """
    if isinstance(norm, mcolors.Normalize):
        return norm
    if isinstance(norm, str):
        # Get class
        norm_out = normalizers.get(norm, None)
        if norm_out is None:
            raise ValueError(f'Unknown normalizer "{norm}". Options are {", ".join(normalizers.keys())}.')
        # Instantiate class
        if norm_out is LinearSegmentedNorm:
            if not np.iterable(levels):
                raise ValueError(f'Need levels for normalizer "{norm}". Received levels={levels}.')
            kwargs.update({'levels':levels})
        norm_out = norm_out(**kwargs) # initialize
    else:
        raise ValueError(f'Unknown norm "{norm_out}".')
    return norm_out

#------------------------------------------------------------------------------
# Meta-normalizer class for discrete levels
#------------------------------------------------------------------------------
# See this post: https://stackoverflow.com/a/48614231/4970632
# WARNING: Many methods in ColorBarBase tests for class membership, crucially
# including _process_values(), which if it doesn't detect BoundaryNorm will
# end up trying to infer boundaries from inverse() method. So make it parent.
class BinNorm(mcolors.BoundaryNorm):
    """
    This normalizer is used for all colormap plots. It can be thought of as a
    "parent" normalizer: it first scales the data according to any
    arbitrary `~matplotlib.colors.Normalize` class, then maps the normalized
    values ranging from 0-1 into **discrete** levels.

    Consider input levels of ``[0, 3, 6, 9, 12, 15]``. The algorithm is
    as follows.

    1. `levels` are normalized according to the input normalizer `norm`.
       If it is ``None``, they are not changed. Possible normalizers include
       `~matplotlib.colors.LogNorm`, which makes color transitions linear in
       the logarithm of the value, or `LinearSegmentedNorm`, which makes
       color transitions linear in the **index** of the level array.
    2. Possible colormap coordinates, corresponding to bins delimited by the
       normalized `levels` array, are calculated.  In this case, the bin centers
       are simply ``[1.5, 4.5, 7.5, 10.5, 13.5]``, which gives us normalized
       colormap coordinates of ``[0, 0.25, 0.5, 0.75, 1]``.
    3. Out-of-bounds coordinates are added. These depend on the value of the
       `extend` keyword argument. For `extend` equal to ``'neither'``,
       the coordinates including out-of-bounds values are ``[0, 0, 0.25, 0.5, 0.75, 1, 1]`` --
       out-of-bounds values have the same color as the nearest in-bounds values.
       For `extend` equal to ``'both'``, the bins are ``[0, 0.16, 0.33, 0.5, 0.66, 0.83, 1]`` --
       out-of-bounds values are given distinct colors. This makes sure your
       colorbar always shows the **full range of colors** in the colormap.
    4. Whenever `BinNorm.__call__` is invoked, the input value normalized by
       `norm` is compared against the normalized `levels` array. Its bin index
       is determined with `numpy.searchsorted`, and its corresponding
       colormap coordinate is selected using this index.

    The input parameters are as follows.
    """
    def __init__(self, levels, norm=None, clip=False, step=1.0, extend='neither'):
        """
        Parameters
        ----------
        levels : list of float
            The discrete data levels.
        norm : `~matplotlib.colors.Normalize`, optional
            The normalizer used to transform `levels` and all data passed
            to `BinNorm.__call__` *before* discretization.
        step : float, optional
            The intensity of the transition to out-of-bounds color, as a
            faction of the *average* step between in-bounds colors.
            Default is ``1``.
        extend : {'neither', 'both', 'min', 'max'}, optional
            Which direction colors will be extended. No matter the `extend`
            option, `BinNorm` ensures colors always extend through the
            extreme end colors.
        clip : bool, optional
            A `~matplotlib.colors.Normalize` option.

        Note
        ----
        If you are using a diverging colormap with ``extend='max'`` or
        ``extend='min'``, the center will get messed up. But that is very strange
        usage anyway... so please just don't do that :)
        """
        # Declare boundaries, vmin, vmax in True coordinates.
        # Notes:
        # * Idea is that we bin data into len(levels) discrete x-coordinates,
        #   and optionally make out-of-bounds colors the same or different
        # * Don't need to call parent __init__, this is own implementation
        #   Do need it to subclass BoundaryNorm, so ColorbarBase will detect it
        #   See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        levels = np.atleast_1d(levels)
        if levels.size <= 1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1]) <= 0).any():
            raise ValueError(f'Levels {levels} passed to Normalize() must be monotonically increasing.')
        if extend not in ('both','min','max','neither'):
            raise ValueError(f'Unknown extend option "{extend}". Choose from "min", "max", "both", "neither".')

        # Determine color ids for levels, i.e. position in 0-1 space
        # Length of these ids should be N + 1 -- that is, N - 1 colors
        # for values in-between levels, plus 2 colors for out-of-bounds.
        # * For same out-of-bounds colors, looks like [0, 0, ..., 1, 1]
        # * For unique out-of-bounds colors, looks like [0, X, ..., 1 - X, 1]
        #   where the offset X equals step/len(levels).
        # First get coordinates
        if not norm:
            norm = mcolors.Normalize() # WARNING: Normalization to 0-1 must always take place first, required by colorbar_factory ticks manager.
        x_b = norm(levels)
        x_m = (x_b[1:] + x_b[:-1])/2 # get level centers after norm scaling
        y = (x_m - x_m.min())/(x_m.max() - x_m.min())
        if isinstance(y, ma.core.MaskedArray):
            y = y.filled(np.nan)
        y = y[np.isfinite(y)]
        # Account for out of bounds colors
        # WARNING: For some reason must clip manually for LogNorm, or
        # end up with unpredictable fill value, weird "out-of-bounds" colors
        offset = 0
        scale = 1
        eps = step/(y.size - 1)
        # eps = step/levels.size
        if extend in ('min','both'):
            offset = eps
            scale -= eps
        if extend in ('max','both'):
            scale -= eps
        y = np.concatenate(([0], offset + scale*y, [1])) # insert '0' (arg 3) before index '0' (arg 2)
        self._norm = norm
        self._x_b = x_b
        self._y = y
        if isinstance(norm, mcolors.LogNorm):
            self._norm_clip = (5e-249, None)
        else:
            self._norm_clip = None

        # Add builtin properties
        # NOTE: Are vmin/vmax even used?
        self.boundaries = levels
        self.vmin = levels.min()
        self.vmax = levels.max()
        self.clip = clip
        self.N = levels.size

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1."""
        # Follow example of LinearSegmentedNorm, but perform no interpolation,
        # just use searchsorted to bin the data.
        norm_clip = self._norm_clip
        if norm_clip:
            xq = np.clip(xq, *norm_clip)
        xq = self._norm(xq)
        yq = self._y[np.searchsorted(self._x_b, xq)] # which x-bin does each point in xq belong to?
        mask = ma.getmaskarray(xq)
        return ma.array(yq, mask=mask)

    def inverse(self, yq):
        """Raises error -- inversion after discretization is impossible."""
        raise ValueError('BinNorm is not invertible.')

#-----------------------------------------------------------------------------#
# Normalizers intended to *pre-scale* levels passed to BinNorm
#-----------------------------------------------------------------------------#
class LinearSegmentedNorm(mcolors.Normalize):
    """
    This is the default normalizer paired with `BinNorm` whenever `levels`
    are non-linearly spaced. The normalized value is linear with respect to
    its **average index** in the `levels` vector, allowing uniform color transitions
    across **arbitrarily spaced** monotonically increasing values.

    It accomplishes this following the example of the `~matplotlib.colors.LinearSegmentedColormap`
    source code, by performing efficient, vectorized linear interpolation
    between the provided boundary levels.

    Can be used by passing ``norm='segmented'`` or ``norm='segments'`` to any
    command accepting ``cmap``. The default midpoint is zero.
    """
    def __init__(self, levels, vmin=None, vmax=None, **kwargs):
        """
        Parameters
        ----------
        levels : list of float
            The discrete data levels.
        vmin, vmax : None
            Ignored, because `vmin` and `vmax` are set to the minimum and
            maximum of `levels`.
        **kwargs
            Passed to `~matplotlib.colors.Normalize`.
        """
        # Save levels
        levels = np.atleast_1d(levels)
        if levels.size <= 1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1]) <= 0).any():
            raise ValueError(f'Levels {levels} passed to LinearSegmentedNorm must be monotonically increasing.')
        vmin, vmax = levels.min(), levels.max()
        super().__init__(vmin, vmax, **kwargs) # second level superclass
        self._x = levels
        self._y = np.linspace(0, 1, len(levels))

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1. Inverse operation
        of `~LinearSegmentedNorm.inverse`."""
        # Follow example of make_mapping_array for efficient, vectorized
        # linear interpolation across multiple segments.
        # Notes:
        # * Normal test puts values at a[i] if a[i-1] < v <= a[i]; for
        #   left-most data, satisfy a[0] <= v <= a[1]
        # * searchsorted gives where xq[i] must be inserted so it is larger
        #   than x[ind[i]-1] but smaller than x[ind[i]]
        x = self._x # from arbitrarily spaced monotonic levels
        y = self._y # to linear range 0-1
        xq = np.atleast_1d(xq)
        ind = np.searchsorted(x, xq)
        ind[ind == 0] = 1
        ind[ind == len(x)] = len(x) - 1 # actually want to go to left of that
        distance = (xq - x[ind - 1])/(x[ind] - x[ind - 1])
        yq = distance*(y[ind] - y[ind - 1]) + y[ind - 1]
        mask = ma.getmaskarray(xq)
        return ma.array(yq, mask=mask)

    def inverse(self, yq):
        """Inverse operation of `~LinearSegmentedNorm.__call__`."""
        x = self._x
        y = self._y
        yq = np.atleast_1d(yq)
        ind = np.searchsorted(y, yq)
        ind[ind == 0] = 1
        ind[ind == len(y)] = len(y) - 1
        distance = (yq - y[ind - 1])/(y[ind] - y[ind - 1])
        xq = distance*(x[ind] - x[ind - 1]) + x[ind - 1]
        mask = ma.getmaskarray(yq)
        return ma.array(xq, mask=mask)

class MidpointNorm(mcolors.Normalize):
    """
    Ensures a "midpoint" always lies at the central colormap color.
    Can be used by passing ``norm='midpoint'`` to any command accepting
    ``cmap``. The default midpoint is zero.
    """
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=None):
        """
        Parameters
        ----------
        midpoint : float, optional
            The midpoint, or the data value corresponding to the normalized
            value ``0.5`` -- halfway down the colormap.
        vmin, vmax, clip
            The minimum and maximum data values, and the clipping setting.
            Passed to `~matplotlib.colors.Normalize`.

        Note
        ----
        See `this stackoverflow thread <https://stackoverflow.com/q/25500541/4970632>`_.
        """
        # Bigger numbers are too one-sided
        super().__init__(vmin, vmax, clip)
        self._midpoint = midpoint

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1. Inverse operation of
        `~MidpointNorm.inverse`."""
        # Get middle point in 0-1 coords, and value
        # Notes:
        # * Look up these three values in case vmin/vmax changed; this is
        #   a more general normalizer than the others. Others are 'parent'
        #   normalizers, meant to be static more or less.
        # * searchsorted gives where xq[i] must be inserted so it is larger
        #   than x[ind[i]-1] but smaller than x[ind[i]]
        #   x, y = [self.vmin, self._midpoint, self.vmax], [0, 0.5, 1]
        if self.vmin >= self._midpoint or self.vmax <= self._midpoint:
            raise ValueError(f'Midpoint {self._midpoint} outside of vmin {self.vmin} and vmax {self.vmax}.')
        x = np.array([self.vmin, self._midpoint, self.vmax])
        y = np.array([0, 0.5, 1])
        xq = np.atleast_1d(xq)
        ind = np.searchsorted(x, xq)
        ind[ind == 0] = 1 # in this case will get normed value <0
        ind[ind == len(x)] = len(x) - 1 # in this case, will get normed value >0
        distance = (xq - x[ind - 1])/(x[ind] - x[ind - 1])
        yq = distance*(y[ind] - y[ind - 1]) + y[ind - 1]
        mask = ma.getmaskarray(xq)
        return ma.array(yq, mask=mask)

    def inverse(self, yq, clip=None):
        """Inverse operation of `~MidpointNorm.__call__`."""
        # Invert the above
        # x, y = [self.vmin, self._midpoint, self.vmax], [0, 0.5, 1]
        # return ma.masked_array(np.interp(yq, y, x))
        # Performs inverse operation of __call__
        x = np.array([self.vmin, self._midpoint, self.vmax])
        y = np.array([0, 0.5, 1])
        yq = np.atleast_1d(yq)
        ind = np.searchsorted(y, yq)
        ind[ind == 0] = 1
        ind[ind == len(y)] = len(y) - 1
        distance = (yq - y[ind - 1])/(y[ind] - y[ind - 1])
        xq = distance*(x[ind] - x[ind - 1]) + x[ind - 1]
        mask = ma.getmaskarray(yq)
        return ma.array(xq, mask=mask)

#-----------------------------------------------------------------------------#
# Register data
#-----------------------------------------------------------------------------#
def _read_cmap_cycle_data(filename):
    """
    Helper function that reads generalized colormap and color cycle files.
    """
    empty = (None, None, None)
    if os.path.isdir(filename): # no warning
        return empty
    # Directly read segmentdata json file
    # NOTE: This is special case! Immediately return name and cmap
    split = os.path.basename(filename).split('.')
    if len(split) == 1:
        return empty
    *name, ext = split
    name = ''.join(name)
    if ext == 'json':
        with open(filename, 'r') as f:
            data = json.load(f)
        N = rcParams['image.lut']
        if 'red' in data:
            cmap = mcolors.LinearSegmentedColormap(name, data, N=N)
        else:
            kw = {}
            for key in ('space','gamma1','gamma2'):
                kw[key] = data.pop(key, None)
            cmap = PerceptuallyUniformColormap(name, data, N=N, **kw)
        if name[-2:] == '_r':
            cmap = cmap.reversed(name[:-2])
        return name, None, cmap
    # Read .rgb, .rgba, .xrgb, and .xrgba files
    elif ext in ('rgb', 'xrgb', 'rgba', 'xrgba'):
        # Load
        # NOTE: This appears to be biggest import time bottleneck! Increases
        # time from 0.05s to 0.2s, with numpy loadtxt or with this regex thing.
        delim = re.compile('[,\s]+')
        data = [delim.split(line.strip()) for line in open(filename).readlines()]
        try:
            data = [[float(num) for num in line] for line in data]
        except ValueError:
            warnings.warn(f'Failed to load "{filename}". Expected a table of comma or space-separated values.')
            return empty
        # Build x-coordinates and standardize shape
        data = np.array(data)
        if data.shape[1] != len(ext):
            warnings.warn(f'Failed to load "{filename}". Got {data.shape[1]} columns, but expected {len(ext)}.')
            return empty
        if ext[0] != 'x': # i.e. no x-coordinates specified explicitly
            x = np.linspace(0, 1, data.shape[0])
        else:
            x, data = data[:,0], data[:,1:]
    # Load XML files created with scivizcolor
    # Adapted from script found here: https://sciviscolor.org/matlab-matplotlib-pv44/
    elif ext == 'xml':
        try:
            xmldoc = etree.parse(filename)
        except IOError:
            warnings.warn(f'Failed to load "{filename}".')
            return empty
        x, data = [], []
        for s in xmldoc.getroot().findall('.//Point'):
            # Verify keys
            if any(key not in s.attrib for key in 'xrgb'):
                warnings.warn(f'Failed to load "{filename}". Missing an x, r, g, or b specification inside one or more <Point> tags.')
                return empty
            if 'o' in s.attrib and 'a' in s.attrib:
                warnings.warn(f'Failed to load "{filename}". Contains ambiguous opacity key.')
                return empty
            # Get data
            color = []
            for key in 'rgbao': # o for opacity
                if key not in s.attrib:
                    continue
                color.append(float(s.attrib[key]))
            x.append(float(s.attrib['x']))
            data.append(color)
        # Convert to array
        if not all(len(data[0]) == len(color) for color in data):
             warnings.warn(f'File {filename} has some points with alpha channel specified, some without.')
             return empty
    elif ext == 'hex':
        # Read hex strings
        string = open(filename).read() # into single string
        data = re.findall('#[0-9a-fA-F]{6}', string) # list of strings
        if len(data) < 2:
            warnings.warn(f'Failed to load "{filename}".')
            return empty
        # Convert to array
        x = np.linspace(0, 1, len(data))
        data = [to_rgb(color) for color in data]
    else:
        warnings.warn(f'Colormap/cycle file "{filename}" has unknown extension.')
        return empty
    # Standardize and reverse if necessary to cmap
    x, data = np.array(x), np.array(data)
    x = (x - x.min()) / (x.max() - x.min()) # for some reason, some aren't in 0-1 range
    if (data > 2).any(): # from 0-255 to 0-1
        data = data/255
    if name[-2:] == '_r':
        name = name[:-2]
        data = data[::-1,:]
        x = 1 - x[::-1]
    # Return data
    return name, x, data

@_timer
def register_cmaps():
    """
    Adds colormaps packaged with ProPlot or saved to the ``~/.proplot/cmaps``
    folder. This is called on import. Maps are registered according to their
    filenames -- for example, ``name.xyz`` will be registered as ``'name'``.
    Use `show_cmaps` to generate a table of the registered colormaps

    Valid file formats are described in the below table.

    =====================  =============================================================================================================================================================================================================
    Extension              Description
    =====================  =============================================================================================================================================================================================================
    ``.hex``               List of HEX strings in any format (comma-separated, separate lines, with double quotes... anything goes).
    ``.xml``               XML files with ``<Point .../>`` entries specifying ``x``, ``r``, ``g``, ``b``, and optionally, ``a`` values, where ``x`` is the colormap coordinate and the rest are the RGB and opacity (or "alpha") values.
    ``.rgb``               3-column table delimited by commas or consecutive spaces, each column indicating red, blue and green color values.
    ``.xrgb``              As with ``.rgb``, but with 4 columns. The first column indicates the colormap coordinate.
    ``.rgba``, ``.xrgba``  As with ``.rgb``, ``.xrgb``, but with a trailing opacity (or "alpha") column.
    =====================  =============================================================================================================================================================================================================
    """
    # Turn original matplotlib maps from ListedColormaps to LinearSegmentedColormaps
    # It makes zero sense to me that they are stored as ListedColormaps
    _check_data()
    for name in CMAPS_CATEGORIES['Matplotlib Originals']: # initialize as empty lists
        cmap = mcm.cmap_d._getitem(name, None)
        if cmap and isinstance(cmap, mcolors.ListedColormap):
            mcm.cmap_d[name] = mcolors.LinearSegmentedColormap.from_list(name, cmap.colors)

    # Misc tasks
    cmap = mcm.cmap_d.pop('Greys', None)
    if cmap is not None:
        mcm.cmap_d['Grays'] = cmap # to be consistent with registered color names (also 'Murica)
    for name in ('Spectral',):
        mcm.cmap_d[name] = mcm.cmap_d[name].reversed(name=name) # make spectral go from 'cold' to 'hot'

    # Remove gross cmaps (strong-arm user into using the better ones)
    for name in CMAPS_DELETE:
        mcm.cmap_d.pop(name, None)

    # Fill initial user-accessible cmap list with the colormaps we will keep
    cmaps.clear()
    cmaps[:] = [name for name,cmap in mcm.cmap_d.items() if not isinstance(cmap,mcolors.ListedColormap)]

    # Add colormaps from ProPlot and user directories
    N = rcParams['image.lut'] # query this when register function is called
    for filename in sorted(glob.glob(os.path.join(DATA_CMAPS, '*'))) + \
            sorted(glob.glob(os.path.join(DATA_USER_CMAPS, '*'))):
        name, x, data = _read_cmap_cycle_data(filename)
        if name is None:
            continue
        if isinstance(data, mcolors.LinearSegmentedColormap):
            cmap = data
        else:
            data = [(x,color) for x,color in zip(x,data)]
            cmap = mcolors.LinearSegmentedColormap.from_list(name, data, N=N)
        mcm.cmap_d[name] = cmap
        cmaps.append(name)
    # Add cyclic attribute
    for name,cmap in mcm.cmap_d.items():
        cmap._cyclic = (name.lower() in ('twilight', 'twilight_shifted', 'phase', 'graycycle')) # add hidden attribute used by BinNorm

    # Sort
    cmaps[:] = sorted(cmaps)

@_timer
def register_cycles():
    """
    Adds color cycles packaged with ProPlot or saved to the ``~/.proplot/cycles``
    folder. This is called on import. Cycles are registered according to their
    filenames -- for example, ``name.hex`` will be registered under the name
    ``'name'`` as a `~matplotlib.colors.ListedColormap` map (see `Cycle` for
    details). Use `show_cycles` to generate a table of the registered cycles.

    For valid file formats, see `register_cmaps`.
    """
    # Empty out user-accessible cycle list
    _check_data()
    cycles.clear()

    # Remove gross cycles, change the names of some others
    for name in CYCLES_DELETE:
        mcm.cmap_d.pop(name, None)
    for (name1,name2) in CYCLES_RENAME:
        cycle = mcm.cmap_d.pop(name1, None)
        if cycle:
            mcm.cmap_d[name2] = cycle
            cycles.append(name2)

    # Read cycles from directories
    icycles = {}
    for filename in sorted(glob.glob(os.path.join(DATA_CYCLES, '*'))) + \
            sorted(glob.glob(os.path.join(DATA_USER_CYCLES, '*'))):
        name, _, data = _read_cmap_cycle_data(filename)
        if name is None:
            continue
        if isinstance(data, mcolors.LinearSegmentedColormap):
            warnings.warn(f'Failed to load {filename} as color cycle.')
            continue
        icycles[name] = data

    # Register cycles as ListedColormaps
    for name,colors in {**CYCLES_PRESET, **icycles}.items():
        cmap = mcolors.ListedColormap(colors, name=name)
        cmap.colors = [to_rgb(color) for color in cmap.colors] # sanitize
        mcm.cmap_d[name] = cmap
        cycles.append(name)

    # Sort
    cycles[:] = sorted([*cycles, 'Set2', 'Set3'], key=lambda s: s.lower())

@_timer
def register_colors(nmax=np.inf):
    """
    Reads full database of crowd-sourced XKCD color names and official
    Crayola color names, then filters them to be sufficiently "perceptually
    distinct" in the HCL colorspace. This is called on import. Use `show_colors`
    to generate a table of the resulting filtered colors.
    """
    # Reset native colors dictionary and add some default groups
    # Add in CSS4 so no surprises for user, but we will not encourage this
    # usage and will omit CSS4 colors from the demo table.
    scale = (360, 100, 100)
    base = {}
    colordict.clear()
    base.update(mcolors.BASE_COLORS)
    base.update(BASE_COLORS_FULL)
    mcolors.colorConverter.colors.clear() # clean out!
    mcolors.colorConverter.cache.clear() # clean out!
    for name,dict_ in (('base',base), ('css',mcolors.CSS4_COLORS)):
        colordict.update({name:dict_})

    # Load colors from file and get their HCL values
    names = ('opencolors', 'xkcd', 'crayola') # order is preference for identical color names from different groups
    files = [os.path.join(DATA_COLORS, f'{name}.txt') for name in names]
    pairs = []
    seen = {*base} # never overwrite base names, e.g. 'blue' and 'b'!
    hcls = np.empty((0,3))
    for file in files:
        category, _ = os.path.splitext(os.path.basename(file))
        data = np.genfromtxt(file, delimiter='\t', dtype=str, comments='%', usecols=(0,1)).tolist()
        # Immediately add all opencolors
        if category == 'opencolors':
            dict_ = {name:color for name,color in data}
            colordict.update({'opencolors':dict_})
            continue
        # Other color dictionaries are filtered, and their names are sanitized
        i = 0
        dict_ = {}
        ihcls = []
        colordict[category] = {} # just initialize this one
        for name,color in data: # is list of name, color tuples
            if i >= nmax: # e.g. for xkcd colors
                break
            for regex,sub in FILTER_TRANSLATIONS:
                name = regex.sub(sub, name)
            if name in seen or FILTER_BAD.search(name):
                continue
            seen.add(name)
            pairs.append((category, name)) # save the category name pair
            ihcls.append(to_xyz(color, space=FILTER_SPACE))
            dict_[name] = color # save the color
            i += 1
        _colordict_unfiltered[category] = dict_
        hcls = np.concatenate((hcls, ihcls), axis=0)

    # Remove colors that are 'too similar' by rounding to the nearest n units
    # WARNING: Unique axis argument requires numpy version >=1.13
    deleted = 0
    hcls = hcls/np.array(scale)
    hcls = np.round(hcls/FILTER_THRESH).astype(np.int64)
    _, idxs, _ = np.unique(hcls, return_index=True, return_counts=True, axis=0) # get unique rows
    for idx,(category,name) in enumerate(pairs):
        if name not in FILTER_ADD and idx not in idxs:
            deleted += 1
        else:
            colordict[category][name] = _colordict_unfiltered[category][name]
    # Add to colors mapping
    for _,kw in colordict.items():
        mcolors.colorConverter.colors.update(kw)

#-----------------------------------------------------------------------------#
# Font stuff
#-----------------------------------------------------------------------------#
# Helvetica: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Valid styles: https://matplotlib.org/api/font_manager_api.html for valid weights, styles, etc.
# Classic fonts: https://www.lifewire.com/classic-sans-serif-fonts-clean-appearance-1077406
# For downloading fonts: https://www.cufonfonts.com
# Notes on getting ttf files on Mac
# * Location in /System/Library/Font, /Library/Fonts, or ~/Library/Fonts
# * To break down .dfont files, use fondu (homebrew download).
#   To break down .ttc files, use dfontsplitter (https://peter.upfold.org.uk/projects/dfontsplitter)
#   To break down .bdf files made by fondu, use mkttf (https://github.com/Tblue/mkttf; requires FontForge and PoTrace)
# * Install new fonts with "brew cask install font- < name-of-font > " after using
#   "brew tap caskroom/fonts" to initialize; appear in ~/Library/Fonts; see https://github.com/Homebrew/homebrew-cask-fonts
# * The .otf files work in addition to .ttf files. You can verify this by
#   looking at plot.fonts_files_os -- it includes .otf files.
# Notes on default files packaged in font directory:
# * Location will be something like:
#   /lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf
# * 'STIX' fonts allow different LaTeX-like math modes e.g. blackboard bold
#   and caligraphy; see: https://matplotlib.org/gallery/text_labels_and_annotations/stix_fonts_demo.html
# * The 'cm'-prefix fonts seem to provide additional mathematical symbols
#   like integrals, and italized math-mode fonts.
# * We also have 'pdfcorefonts' in this directory, but I think since these
#   are afm matplotlib cannot use them? Don't know.
# WARNING: Check out ttflist whenever adding new ttf files! For example, realized
# could dump all of the Gotham-Name.ttf files instead of GothamName files, and
# got Helvetica bug due to unrecognized 'thin' font style overwriting normal one.
# print(*[font for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
# print(*[font.fname for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
@_timer
def register_fonts():
    """Adds fonts packaged with ProPlot or saved to the ``~/.proplot/fonts``
    folder. Also deletes the font cache, which may cause delays.
    Detects ``.ttf`` and ``.otf`` files -- see `this link
    <https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting various other font file types to ``.ttf`` and
    ``.otf`` for use with matplotlib."""
    # Add proplot path to TTFLIST and rebuild cache
    # NOTE: Delay font_manager import, because want to avoid rebuilding font
    # cache, which means import must come after TTFPATH added to environ!
    _check_data()
    paths = DATA_FONTS + ':' + DATA_USER_FONTS
    if 'TTFPATH' not in os.environ:
        os.environ['TTFPATH'] = paths
    elif paths not in os.environ['TTFPATH']:
        os.environ['TTFPATH'] += (':' + paths)

    # Load font manager and rebuild only if necessary!
    # Font cache rebuild can be >50% of total import time, ~1s!!!
    import matplotlib.font_manager as mfonts
    files_loaded = {font.fname for font in mfonts.fontManager.ttflist}
    files_ttfpath = {*mfonts.findSystemFonts(paths.split(':'))}
    if not (files_ttfpath <= files_loaded):
        mfonts._rebuild()

    # Populate font lists
    fonts_system[:] = sorted({font.name for font in mfonts.fontManager.ttflist if not (DATA_USER_FONTS in font.fname or DATA_FONTS in font.fname)})
    fonts_proplot[:] =  sorted({font.name for font in mfonts.fontManager.ttflist if (DATA_USER_FONTS in font.fname or DATA_FONTS in font.fname)})
    fonts[:] = sorted((*fonts_system, *fonts_proplot))

#-----------------------------------------------------------------------------#
# Register stuff and define variables
#-----------------------------------------------------------------------------#
cmaps = [] # track *downloaded* colormaps, user can then check this list
"""List of new registered colormap names."""
cycles = [] # track *all* color cycles
"""List of registered color cycle names."""
_colordict_unfiltered = {} # downloaded colors categorized by filename
colordict = {} # limit to 'sufficiently unique' color names
"""Registered color names by category."""
fonts_proplot = []
"""Names of fonts added by ProPlot."""
fonts_system = []
"""Names of fonts provided by matplotlib or your operating system."""
fonts = []
"""All registered font names."""

register_colors()
register_cmaps()
register_cycles()
register_fonts()

# Dictionary of normalizers; note BinNorm is inaccessible for users
normalizers = {
    'none':       mcolors.NoNorm,
    'null':       mcolors.NoNorm,
    'zero':       MidpointNorm,
    'midpoint':   MidpointNorm,
    'segments':   LinearSegmentedNorm,
    'segmented':  LinearSegmentedNorm,
    'log':        mcolors.LogNorm,
    'linear':     mcolors.Normalize,
    'power':      mcolors.PowerNorm,
    'symlog':     mcolors.SymLogNorm,
    }
"""Dictionary of possible normalizers. See `Norm` for a table."""

#-----------------------------------------------------------------------------#
# Demos
#-----------------------------------------------------------------------------#
def show_channels(*args, N=100, rgb=True, minhue=0, width=100,
    aspect=1, axwidth=1.7):
    """
    Shows how arbitrary colormap(s) vary with respect to the hue, chroma,
    luminance, HSL saturation, and HPL saturation channels, and optionally
    the red, blue and green channels. Adapted from
    `this example <https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html#lightness-of-matplotlib-colormaps>`__.

    Parameters
    ----------
    *args : colormap-spec, optional
        Positional arguments are colormap names or objects. Default is
        :rc:`image.cmap`.
    N : int, optional
        The number of markers to draw for each colormap.
    rgb : bool, optional
        Whether to also show the red, blue, and green channels in the bottom
        row. Default is ``True``.
    minhue : float, optional
        The minimum hue. This lets you rotate the hue plot cyclically.
    width : int, optional
        The width of each colormap line in points.
    aspect : float or (float,float), optional
        The aspect ratio of the subplot.
    axwidth : float, optional
        The width of the subplots.

    Returns
    -------
    `~proplot.subplots.Figure`
        The figure object.
    """
    # Figure and plot
    from . import subplots
    if not args:
        args = (rcParams['image.cmap'],)
    array = [[1,1,2,2,3,3],[0,4,4,5,5,0],[6,6,7,7,8,8]]
    if not rgb:
        array = array[:2]
    fig, axs = subplots(
        array=array, axwidth=axwidth, span=False, share=1,
        aspect=aspect, axpad='1em',
        )
    labels = (
        'Hue', 'Chroma', 'Luminance',
        'HSL saturation', 'HPL saturation',
        'Red', 'Blue', 'Green'
        )
    # Iterate through colormaps
    mc, ms, mp = 0, 0, 0
    cmaps = []
    for cmap in args:
        # Get colormap
        cmap = Colormap(cmap, N=N) # arbitrary cmap argument
        cmap._init()
        cmaps.append(cmap)
        # Get clipped RGB table
        x = np.linspace(0, 1, N)
        lut = cmap._lut[:-3,:3].copy()
        rgb = lut.T # 3 by N
        hcl = np.array([to_xyz(color, space='hcl') for color in lut]).T # 3 by N
        hsl = [to_xyz(color, space='hsl')[1] for color in lut]
        hpl = [to_xyz(color, space='hpl')[1] for color in lut]
        # Plot channels
        # If rgb is False, the zip will just truncate the other iterables
        for ax,y,label in zip(axs,(*hcl,hsl,hpl,*rgb),labels):
            ylim, ylocator = None, None
            if label in ('Red','Blue','Green'):
                ylim = (0,1)
                ylocator = 0.2
            elif label == 'Luminance':
                ylim = (0,100)
                ylocator = 20
            elif label == 'Hue':
                ylim = (minhue,minhue+360)
                ylocator = 90
                y = y - 720
                for i in range(3): # rotate up to 1080 degrees
                    y[y < minhue] += 360
            else:
                if label == 'Chroma':
                    mc = max(mc,max(y))
                    m = mc
                elif 'HSL' in label:
                    ms = max(ms,max(y))
                    m = ms
                else:
                    mp = max(mp,max(y))
                    m = mp
                ylim = (0,m)
                ylocator = ('maxn', 5)
            ax.scatter(x, y, c=x, cmap=cmap, s=width, linewidths=0)
            ax.format(title=label, ylim=ylim, ylocator=ylocator)
    # Formatting
    suptitle = ', '.join(f'"{cmap.name}"' for cmap in cmaps[:-1]) \
        + (', and ' if len(cmaps) > 2 else ' and ' if len(cmaps) == 2 else ' ') \
        + f'"{cmaps[-1].name}" colormap' + ('s' if len(cmaps) > 1 else '')
    axs.format(
        xlocator=0.25, xformatter='null',
        suptitle=f'{suptitle} by channel', ylim=None, ytickminor=False,
        )
    # Colorbar on the bottom
    for cmap in cmaps:
        fig.colorbar(cmap,
            loc='b', span=(2,5),
            locator='null', label=cmap.name, labelweight='bold')
    return fig

def show_colorspaces(luminance=None, saturation=None, hue=None):
    """
    Generates hue-saturation, hue-luminance, and luminance-saturation
    cross-sections for the HCL, HSLuv, and HPLuv colorspaces.

    Parameters
    ----------
    luminance : float, optional
        If passed, chroma-saturation cross-sections are drawn for this luminance.
        Must be between ``0` and ``100``. Default is ``50``.
    saturation : float, optional
        If passed, luminance-hue cross-sections are drawn for this saturation.
        Must be between ``0` and ``100``.
    hue : float, optional
        If passed, luminance-saturation cross-sections are drawn for this hue.
        Must be between ``0` and ``360``.

    Returns
    -------
    `~proplot.subplots.Figure`
        The figure object.
    """
    # Get colorspace properties
    hues = np.linspace(0, 360, 361)
    sats = np.linspace(0, 120, 120) # use 120 instead of 121, prevents annoying rough edge on HSL plot
    lums = np.linspace(0, 99.99, 101)
    if luminance is None and saturation is None and hue is None:
        luminance = 50
    if luminance is not None:
        hsl = np.concatenate((
            np.repeat(hues[:,None], len(sats), axis=1)[...,None],
            np.repeat(sats[None,:], len(hues), axis=0)[...,None],
            np.ones((len(hues), len(sats)))[...,None]*luminance,
            ), axis=2)
        suptitle = f'Hue-saturation cross-section for luminance {luminance}'
        xlabel, ylabel = 'hue', 'saturation'
        xloc, yloc = 60, 20
    elif saturation is not None:
        hsl = np.concatenate((
            np.repeat(hues[:,None], len(lums), axis=1)[...,None],
            np.ones((len(hues), len(lums)))[...,None]*saturation,
            np.repeat(lums[None,:], len(hues), axis=0)[...,None],
            ), axis=2)
        suptitle = f'Hue-luminance cross-section for saturation {saturation}'
        xlabel, ylabel = 'hue', 'luminance'
        xloc, yloc = 60, 20
    elif hue is not None:
        hsl = np.concatenate((
            np.ones((len(lums), len(sats)))[...,None]*hue,
            np.repeat(sats[None,:], len(lums), axis=0)[...,None],
            np.repeat(lums[:,None], len(sats), axis=1)[...,None],
            ), axis=2)
        suptitle = 'Luminance-saturation cross-section'
        xlabel, ylabel = 'luminance', 'saturation'
        xloc, yloc = 20, 20

    # Make figure, with black indicating invalid values
    # Note we invert the x-y ordering for imshow
    from . import subplots
    fig, axs = subplots(
        ncols=3, share=0, axwidth=2, aspect=1, axpad=0.05
        )
    for ax,space in zip(axs,('hcl','hsl','hpl')):
        rgba = np.ones((*hsl.shape[:2][::-1], 4)) # RGBA
        for j in range(hsl.shape[0]):
            for k in range(hsl.shape[1]):
                rgb_jk = to_rgb(hsl[j,k,:].flat, space)
                if not all(0 <= c <= 1 for c in rgb_jk):
                    rgba[k,j,3] = 0 # black cell
                else:
                    rgba[k,j,:3] = rgb_jk
        ax.imshow(rgba, origin='lower', aspect='auto')
        ax.format(xlabel=xlabel, ylabel=ylabel, suptitle=suptitle,
                  grid=False, xtickminor=False, ytickminor=False,
                  xlocator=xloc, ylocator=yloc, facecolor='k',
                  title=space.upper(), titleweight='bold')
    return fig

def show_colors(nbreak=17, minsat=0.2):
    """
    Visualizes the registered color names in two figures. Adapted from
    `this example <https://matplotlib.org/examples/color/named_colors.html>`_.

    Parameters
    ----------
    nbreak : int, optional
        The number of breaks between hues for grouping "like colors" in the
        color table.
    minsat : float, optional
        The threshold saturation, between ``0`` and ``1``, for designating
        "gray colors" in the color table.

    Returns
    -------
    (`~proplot.subplots.Figure`, `~proplot.subplots.Figure`)
        A list of the figure objects.
    """
    # Get colors explicitly defined in colorConverter, or the default
    # components of that map
    figs = []
    from . import subplots
    for opencolors in (True,False):
        scale = (360, 100, 100)
        if opencolors:
            group = ['opencolors']
        else:
            group = [name for name in colordict if name not in ('css','opencolors')]
        icolors = {}
        for name in group:
            icolors.update(colordict[name]) # add category dictionary

        # Group colors together by discrete range of hue, then sort by value
        # For opencolors this is not necessary
        if opencolors:
            wscale = 0.5
            swatch = 1.5
            nrows, ncols = 10, len(OPEN_COLORS) # rows and columns
            plot_names = [[name + str(i) for i in range(nrows)] for name in OPEN_COLORS]
            nrows = nrows*2
            ncols = (ncols+1)//2
            plot_names = np.array(plot_names, order='C')
            plot_names.resize((ncols, nrows))
            plot_names = plot_names.tolist()
        # Get colors in perceptally uniform space, then group based on hue thresholds
        else:
            # Transform to HCL space
            ncols = 4
            wscale = 1
            swatch = 1
            colors_hcl = {
                key: [c/s for c,s in zip(to_xyz(value, FILTER_SPACE), scale)]
                for key,value in icolors.items()
                }
            # Separate into columns and roughly sort by brightness in these columns
            breakpoints = np.linspace(0,1,nbreak) # group in blocks of 20 hues
            plot_names = [] # initialize
            sat_test = (lambda x: x < minsat) # test saturation for 'grays'
            for n in range(nbreak):
                # 'Grays' column
                if n == 0:
                    hue_colors = [(name,hcl) for name,hcl in colors_hcl.items() if sat_test(hcl[1])]
                # Column for nth color
                else:
                    b1, b2 = breakpoints[n-1], breakpoints[n]
                    hue_test = ((lambda x: b1 <= x <= b2) if b2
                        is breakpoints[-1] else (lambda x: b1 <= x < b2))
                    hue_colors = [(name,hcl) for name,hcl
                        in colors_hcl.items() if hue_test(hcl[0])
                        and not sat_test(hcl[1])] # grays have separate category
                # Get indices to build sorted list, then append sorted list
                sorted_index = np.argsort([pair[1][2] for pair in hue_colors])
                plot_names.append([hue_colors[i][0] for i in sorted_index])
            # Concatenate those columns so get nice rectangle
            names = [i for sublist in plot_names for i in sublist]
            plot_names = [[]]
            nrows = len(names)//ncols+1
            for i,name in enumerate(names):
                if ((i + 1) % nrows) == 0:
                    plot_names.append([]) # add new empty list
                plot_names[-1].append(name)

        # Create plot by iterating over columns
        fig, ax = subplots(
            width=8*wscale*(ncols/4), height=5*(nrows/40),
            left=0, right=0, top=0, bottom=0, tight=False
            )
        X, Y = fig.get_dpi()*fig.get_size_inches() # size in *dots*; make these axes units
        hsep, wsep = Y/(nrows+1), X/ncols # height and width of row/column in *dots*
        for col,huelist in enumerate(plot_names):
            for row,name in enumerate(huelist): # list of colors in hue category
                if not name: # empty slot
                    continue
                y = Y - hsep*(row + 1)
                y_line = y + hsep*0.1
                xi_line = wsep*(col + 0.05)
                xf_line = wsep*(col + 0.25*swatch)
                xi_text = wsep*(col + 0.25*swatch + 0.03*swatch)
                ax.text(xi_text, y, re.sub('^xkcd:', '', name),
                        fontsize=hsep*0.8, ha='left', va='center')
                ax.hlines(y_line, xi_line, xf_line, color=icolors[name], lw=hsep*0.6)
        # Apply formatting
        ax.format(xlim=(0,X), ylim=(0,Y))
        ax.set_axis_off()
        figs.append(fig)
    return figs

def show_cmaps(*args, N=256, length=4.0, width=0.2):
    """
    Visualizes all registered colormaps, or the list of colormap names if
    positional arguments are passed. Adapted from `this example
    <http://matplotlib.org/examples/color/colormaps_reference.html>`__.

    Parameters
    ----------
    *args : colormap-spec, optional
        Positional arguments are colormap names or objects. Default is
        all of the registered colormaps.
    N : int, optional
        The number of levels in each colorbar.
    length : float or str, optional
        The length of each colorbar. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.
    width : float or str, optional
        The width of each colorbar. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

    Returns
    -------
    `~proplot.subplots.Figure`
        The figure object.
    """
    # Have colormaps separated into categories
    if args:
        imaps = [Colormap(cmap, N=N).name for cmap in args]
    else:
        imaps = [
            name for name in mcm.cmap_d.keys() if name not in ('vega', 'greys', 'no_name')
            and isinstance(mcm.cmap_d[name], mcolors.LinearSegmentedColormap)
            ]

    # Get dictionary of registered colormaps and their categories
    imaps = [name.lower() for name in imaps]
    cats = {cat:names for cat,names in CMAPS_CATEGORIES.items()}
    cats_plot = {cat:[name for name in names if name.lower() in imaps] for cat,names in cats.items()}
    # Distinguish known from unknown (i.e. user) maps, add as a new category
    imaps_known = [name.lower() for cat,names in cats.items() for name in names if name.lower() in imaps]
    imaps_user = [name for name in imaps if name not in imaps_known]
    cats_plot['User'] = imaps_user
    # Remove categories with no known maps
    cats_plot = {cat:maps for cat,maps in cats_plot.items() if maps}

    # Array for producing visualization with imshow
    a = np.linspace(0, 1, 257).reshape(1,-1)
    a = np.vstack((a,a))
    # Figure
    from . import subplots
    naxs = len(imaps_known) + len(imaps_user) + len(cats_plot)
    fig, axs = subplots(
        nrows=naxs, axwidth=length, axheight=width,
        share=0, hspace=0.03,
        )
    iax = -1
    ntitles, nplots = 0, 0 # for deciding which axes to plot in
    for cat,names in cats_plot.items():
        # Space for title
        if not names:
            continue
        ntitles += 1
        for imap,name in enumerate(names):
            # Draw colorbar
            iax += 1
            if imap + ntitles + nplots > naxs:
                break
            ax = axs[iax]
            if imap == 0: # allocate this axes for title
                iax += 1
                ax.set_visible(False)
                ax = axs[iax]
            if name not in mcm.cmap_d or name.lower() not in imaps: # i.e. the expected builtin colormap is missing
                ax.set_visible(False) # empty space
                continue
            ax.imshow(a, cmap=name, origin='lower', aspect='auto', levels=N)
            ax.format(ylabel=name,
                      ylabel_kw={'rotation':0, 'ha':'right', 'va':'center'},
                      xticks='none',  yticks='none', # no ticks
                      xloc='neither', yloc='neither', # no spines
                      title=(cat if imap == 0 else None))
        # Space for plots
        nplots += len(names)
    return fig

def show_cycles(*args, axwidth=1.5):
    """
    Visualizes all registered color cycles, or the list of cycle names if
    positional arguments are passed.

    Parameters
    ----------
    *args : colormap-spec, optional
        Positional arguments are cycle names or objects. Default is
        all of the registered colormaps.
    axwidth : str or float, optional
        Average width of each subplot. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

    Returns
    -------
    `~proplot.subplots.Figure`
        The figure object.
    """
    # Get the list of cycles
    if args:
        icycles = [colors(cycle) for cycle in args]
    else:
        icycles = {key:mcm.cmap_d[key].colors for key in cycles} # use global cycles variable
    nrows = len(icycles)//3 + len(icycles)%3

    # Create plot
    from . import subplots
    state = np.random.RandomState(528)
    fig, axs = subplots(
        ncols=3, nrows=nrows, aspect=1, axwidth=axwidth,
        sharey=False, sharex=False, axpad=0.05
        )
    for i,(ax,(key,cycle)) in enumerate(zip(axs, icycles.items())):
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
    axs[i+1:].set_visible(False)
    return fig

def show_fonts(fonts=None, size=12):
    """Displays table of the fonts installed by ProPlot or in the user-supplied
    `fonts` list. Use `size` to change the fontsize for fonts shown in the figure."""
    # letters = 'Aa Bb Cc Dd Ee Ff Gg Hh Ii Jj Kk Ll Mm Nn Oo Pp Qq Rr Ss Tt Uu Vv Ww Xx Yy Zz'
    from . import subplots
    fonts = ('DejaVu Sans', *fonts_proplot)
    math = r'(0) + {1} - [2] * <3> / 4,0 $\geq\gg$ 5.0 $\leq\ll$ ~6 $\times$ 7 $\equiv$ 8 $\approx$ 9 $\propto$'
    greek = r'$\alpha\beta$ $\Gamma\gamma$ $\Delta\delta$ $\epsilon\zeta\eta$ $\Theta\theta$ $\kappa\mu\nu$ $\Lambda\lambda$ $\Pi\pi$ $\xi\rho\tau\chi$ $\Sigma\sigma$ $\Phi\phi$ $\Psi\psi$ $\Omega\omega$ !?&#%'
    letters = 'the quick brown fox jumps over a lazy dog\nTHE QUICK BROWN FOX JUMPS OVER A LAZY DOG'
    for weight in ('normal',):
        f, axs = subplots(ncols=1, nrows=len(fonts), space=0, axwidth=4.5, axheight=5.5*size/72)
        axs.format(xloc='neither', yloc='neither', xlocator='null', ylocator='null', alpha=0)
        axs[0].format(title='Fonts demo', titlefontsize=size, titleloc='l', titleweight='bold')
        for i,ax in enumerate(axs):
            font = fonts[i]
            ax.text(0, 0.5, f'{font}: {letters}\n{math}\n{greek}', fontfamily=font,
                    fontsize=size, weight=weight, ha='left', va='center')
    return f

#-----------------------------------------------------------------------------#
# Deleted colormaps and colormap categories
#-----------------------------------------------------------------------------#
# TODO: add examples of how to reconstruct e.g. 'tab20c' on-the-fly
# TODO: add examples of how to reconstruct Wave, Insert, Highlight, and Outlier
# Fabio Crameri
# See: http://www.fabiocrameri.ch/colourmaps.php
# 'Fabio Crameri Sequential': [
#     'Acton', 'Buda', 'Lajolla',
#     'Bamako', 'Nuuk', 'Davos', 'Oslo', 'Devon', 'Tokyo',
#     'Batlow', 'Turku', 'Bilbao', 'Lapaz',
#     ],
# 'Fabio Crameri Diverging': [
#     'Roma', 'Broc', 'Cork',  'Vik', 'Oleron',
#     ],
# Kenneth Moreland
# See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/index.html
# Soft coolwarm from: https://www.kennethmoreland.com/color-advice/
# 'Kenneth Moreland': [
#     'CoolWarm', 'MutedCoolWarm', 'SoftCoolWarm',
#     'BlueTan', 'PurpleOrange', 'CyanMauve', 'BlueYellow', 'GreenRed',
#     ],
# 'Kenneth Moreland Sequential': [
#     'BlackBody', 'Kindlmann', 'ExtendedKindlmann',
#     ],
# Elevation and bathymetry
# 'Geographic': [
#     'Bath1', # from XKCD; see: http://soliton.vm.bytemark.co.uk/pub/cpt-city/xkcd/tn/xkcd-bath.png.index.html
#     'Bath2', # from Tom Patterson; see: http://soliton.vm.bytemark.co.uk/pub/cpt-city/tp/index.html
#     'Bath3', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/ibcso/tn/ibcso-bath.png.index.html
#     'Bath4', # ^^ same
#     'Geography4-1', # mostly ocean
#     'Geography5-4', # range must be -4000 to 5000
#     'Geography1', # from ???
#     'Geography2', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/ngdc/tn/ETOPO1.png.index.html
#     'Geography3', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/mby/tn/mby.png.index.html
#     ],
