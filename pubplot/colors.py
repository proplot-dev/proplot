#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Imports
# WARNING: Can't make this module depend on base.subplots, or end up
# with circular imports! Make plots the old-fashioned way.
#------------------------------------------------------------------------------#
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from matplotlib import rcParams
from cycler import cycler

#------------------------------------------------------------------------------
# Colormap stuff
#------------------------------------------------------------------------------
def cmap(name, samples=None, levels=None, norm=None):
    """
    Just wrapper around get_cmap, makes it accessible without needing pyplot import.
    WARNING: Below fails! Lookup table refers to integer values indexing the
    integer list of colors.
    """
    # Previous incorrect attempt
    # lut = levels
    # if utils.isiterable(levels):
    #     lut = min(levels) + [l/(max(levels)-min(levels)) for l in levels]
    # cmap = mcm.get_cmap(name=name)
    # Hacky way that works
    if norm is None:
        norm = mcolors.Normalize(0,1)
    elif not hasattr(norm,'__call__'):
        raise ValueError('Norm has to be callable.')
    elif not isinstance(norm,mcolors.Normalize):
        raise ValueError('Norm has to be an mcolors.Normalize class or subclass.')
    if samples is not None: # use these as the *centers*!
        samples = np.array(samples)
        samples = norm(samples)
        if samples[1]<samples[0]:
            samples = samples[::-1]
            name = name[:-2] if name.endswith('_r') else f'{name}_r'
        idelta = samples[1]-samples[0]
        fdelta = samples[-1]-samples[-2]
        levels = [samples[0]-idelta/2, *((samples[1:]+samples[:-1])/2), samples[-1]+fdelta/2]
        levels = norm.inverse(levels)
    elif levels[1]<levels[0]: # reversed
        levels = levels[::-1]
        name = name[:-2] if name.endswith('_r') else f'{name}_r'
    m = plt.contourf([0,0], [0,0], np.nan*np.ones((2,2)), cmap=name, levels=levels)
    return m

def cmapfactory(colors, levels=None, extend='neither'):
    """
    ***Inverse of cmapcolors.***
    Generate colormap instance from list of levels and colors.
    * Generally don't want these kinds of colorbars to 'extend', but you can.
    * Object will assume one color between individual levels; for example,
      if levels are [0,0.5,0.6,1], the first and last color-zones will be way thicker.
    * If you want unevenly space intervals but evenly spaced boundaries, use custom
      Norm instead of default.
    * This is *one way* to create a colormap from list of colors; another way is
      to just pass a list of colors to any _cformat method. That method will
      also generate levels that are always equally spaced.
    """
    # Use builtin method
    if levels is None:
        levels = np.linspace(0,1,len(colors)+1)
    if len(levels)!=len(colors)+1:
        raise ValueError(f"Have {len(levels):d} levels and {len(colors):d} colors. Need ncolors+1==nlevels.")
    if extend=='min':
        colors = ['w', *colors] # add dummy color
    elif extend=='max':
        colors = [*colors, 'w'] # add dummy color
    elif extend=='both':
        colors = ['w', *colors, 'w']
    elif extend!='neither':
        raise ValueError("Unknown extend option \"{extend}\".")
    cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend=extend)
    return cmap

def cmapcolors(name, N=None, interp=False, vmin=None, vmax=None, sample='center'):
    """
    ***Inverse of cmapfactory.***
    Get individual colors from a discrete colormap, use for e.g. consecutive lines.
    Can also pass list of colors directly to _cformat, to build a colormap; but this
    is sometimes backwards, because will get colors with cmapcolors then build the
    colormap back up again; kind of silly. Consider changing.
    * First argument is the colormap name, second is the number of colors needed.
    * Optionally can input a list, and the colors will be scaled by their relative
      position in vmin/vmax. This is still useful sometimes.
    * Sample can be 'center', 'left', or 'full'.
    """
    cmap = mcm.get_cmap(name) # the cmap object itself
    nsample = N 
    if N is None or not interp: nsample = cmap.N # use the builtin N
    try: iter(nsample)
    except TypeError:
        if sample=='center':
            samples = np.linspace(1/(2*nsample), 1-1/(2*nsample), nsample) # N samples, from centers
        elif sample=='left':
            samples = np.linspace(0, 1-1/nsample, nsample)
        elif sample=='full':
            samples = np.linspace(0, 1, nsample) # N samples, from edge to edge
        else:
            raise ValueError(f'Unknown sampling setting {sample}. Choose from "center", "left", and "full".')
        colors = [cmap(s) for s in samples] # choose these samples
    else:
        if vmin is None or vmax is None:
            raise ValueError("If you input a vector, you must specify the color range for that data with \"vmin\" and \"vmax\".")
        colors = [cmap((v-vmin)/(vmax-vmin)) for v in nsample]
    colors = colors[:N] # trim in case we had extra points
    return colors

def cmapshow(N=11, ignore=['Cycles','Miscellaneous','Sequential2','Diverging2']):
    """
    Plot all current colormaps, along with their catgories.
    This example comes from the Cookbook on www.scipy.org. According to the
    history, Andrew Straw did the conversion from an old page, but it is
    unclear who the original author is.
    See: http://matplotlib.org/examples/color/colormaps_reference.html
    """
    # Have colormaps separated into categories:
    cmaps = [m for m in plt.colormaps() if not m.endswith('_r')]
    cmaps_known = []
    cmaps_custom = ['ColorPicker','HCL','ColorBrewer','Dave','NCL']
    categories = { **{category:[] for category in cmaps_custom}, # initialize as empty lists
        'Cycles': mcolors.CYCLES_CMAPS, # colormaps that are sequences of colors, not gradations
        'Perceptually Uniform Sequential': ['viridis', 'cividis', 'plasma', 'inferno', 'magma'],
        'Diverging1': ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral'],
        'Diverging2': ['coolwarm', 'bwr', 'seismic'],
        'Sequential1': ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
        'Sequential2': ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper'],
        'Miscellaneous': ['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv', 'spectral',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']}
    # print(categories) # for testing
    for k,v in categories.items(): # big list of all colormaps
        cmaps_known += v # add to this existing list
    for i,cm in enumerate(cmaps): # add to 'Custom' if not in above dictionary
        if cm not in cmaps_known:
            if 'Vega' in cm:
                continue # deprecated
            elif cm[:2]=='hc':
                categories['HCL'].append(cm)
            elif cm[:2]=='nc':
                categories['NCL'].append(cm)
            elif cm[:2]=='cb':
                categories['ColorBrewer'].append(cm)
            elif cm[:2]=='cp':
                categories['ColorPicker'].append(cm)
            elif cm[:4]=='dave':
                categories['Dave'].append(cm)
            else:
                print(f"Colormap \"{cm}\" has unknown category.")
                # categories['Custom'].append(cm)
    # Array for producing visualization with imshow
    a = np.linspace(0, 1, 257).reshape(1,-1)
    a = np.vstack((a,a))
    # Figure
    twidth = 1 # number of axes-widths to allocate for titles
    nmaps = len([item for group in categories.values() for item in group]) + len(categories)*twidth
    # nmaps = len(cmaps) + len(categories)*twidth # title for each category
    for cat in ignore:
        nmaps -= (twidth + len(categories[cat])) # reduce size
    fig = plt.figure(figsize=(5,.3*nmaps))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.15, right=0.99)
    # Make plot
    ntitles, nplots = 0, 0 # for deciding which axes to plot in
    cats = [cat for cat in categories if cat not in ignore]
    for cat in cats:
        # Space for title
        ntitles += twidth # two axes-widths
        for i,m in enumerate(categories[cat]):
            if i+ntitles+nplots>nmaps:
                break
            # Get object
            # cmap = mcm.get_cmap(m, lut=N) # interpolate lookup table by calling m(values)
            try:
                cmap = mcm.get_cmap(m) # use default number of colors
            except ValueError:
                print(f'Warning: colormap {m} not found.')
                continue
            # Draw, and make invisible
            ax = plt.subplot(nmaps,1,i+ntitles+nplots)
            for s in ax.spines.values():
                s.set_visible(False)
            # ax.patch.set_alpha(0)
            ax.set_xticks([])
            ax.set_yticks([])
            # Draw colormap
            ax.imshow(a, aspect='auto', cmap=cmap, origin='lower')
            if cat not in cmaps_custom:
                cmaps_known.remove(m)
            # Category title
            if i==0:
                t = ax.title
                t.set_text(cat) # category name
                t.set_visible(True)
            # Label
            yl = ax.yaxis.label
            yl.set_text(m) # map name
            yl.set_visible(True)
            yl.set_rotation(0)
            yl.set_ha('right')
            yl.set_va('center')
        # Space for plots
        nplots += len(categories[cat])
    fig.subplots_adjust(top=0.99, bottom=0.00, left=0.15, right=0.99)
    # Check
    for cat in ignore:
        for m in categories[cat]:
            if cat not in cmaps_custom:
                cmaps_known.remove(m)
    if len(cmaps_known)>0:
        print(f'Colormaps in dictionary, but not found: {", ".join(cmaps_known)}')
    # Save
    # fig.patch.set_alpha(0)
    filename = f'{os.path.dirname(__file__)}/colormaps.pdf'
    print(f"Saving figure to: {filename}.")
    fig.savefig(filename, bbox_inches='tight')
    return fig

#------------------------------------------------------------------------------
# Normalization classes for mapping data to colors (i.e. colormaps)
# WARNING: Many methods in ColorBarBase tests for class membership, crucially
# including _process_values(), which if it doesn't detect BoundaryNorm will
# end up trying to infer boundaries from inverse() method
#------------------------------------------------------------------------------
class ContinuousNorm(mcolors.Normalize):
    """
    As in BoundaryNorm case, but instead we linearly *interpolate* colors
    between the provided boundary indices. Use this e.g. if you want to
    have the gradation from 0-1 'zoomed in' and gradation from 0-10 'zoomed out'.
    In this case, can control number of colors by setting the *lookup table*
    when declaring colormap.
    """
    def __init__(self, levels=None, midpoint=None, clip=False, ncolors=None, **kwargs):
        # Very simple
        levels = np.atleast_1d(levels)
        if np.any((levels[1:]-levels[:-1])<=0):
            raise ValueError(f'Levels passed to ContinuousNorm must be monotonically increasing.')
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

class DiscreteNorm(mcolors.BoundaryNorm):
    """
    Simple normalizer that *interpolates* from an RGB array at point
    (level_idx/num_levels) along the array, instead of choosing color
    from (transform(level_value)-transform(vmin))/(transform(vmax)-transform(vmin))
    where transform can be linear, logarithmic, etc.
    TODO:
    * Allow this to accept transforms too, which will help prevent level edges
      from being skewed toward left or right in case of logarithmic/exponential data.
    Example: Your levels edges are weirdly spaced [-1000, 100, 0, 100, 1000] or
    even [0, 10, 12, 20, 22], but center "colors" are always at colormap
    coordinates [.2, .4, .6, .8] no matter the spacing; levels just must be monotonic.
    """
    def __init__(self, levels=None, midpoint=None, bins=False, clip=False, ncolors=None, **kwargs):
        # Don't need to call partent initializer, this is own implementation; just need
        # it to be subclass so ColorbarBase methods will detect it
        # See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        levels = np.atleast_1d(levels)
        if np.any((levels[1:]-levels[:-1])<=0):
            raise ValueError('Levels passed to Normalize() must be monotonically increasing.')
        self.N = len(levels)
        self.levels = levels
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
            idxs = np.linspace(0, 1, self.levels.size)
            locs = np.where(v>=self.levels)[0]
            if locs.size==0:
                norm[i] = 0
            elif locs.size==self.levels.size:
                norm[i] = 1
            else:
                cutoff  = locs[-1]
                norm[i] = (idxs[cutoff] + idxs[cutoff+1])/2
        return ma.masked_array(norm, np.isnan(value))

    def inverse(self, norm):
        # Not possible
        raise ValueError('Normalize is not invertible in "bins" mode.')

class StretchNorm(mcolors.Normalize):
    """
    Class that can 'stretch' and 'compress' either side of a colormap about
    some midpoint; proceeds exponentially (exp>0) or logarithmically (exp<0)
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
                    block_width*warp((v - midpoint_scaled)/block_width, self.exp)
                    )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled - 
                        block_width*warp((midpoint_scaled - v)/block_width, self.exp)
                        )
        if self.extend=='both' or self.extend=='max':
            value_cmap[value_cmap>1] = 1
        if self.extend=='both' or self.extend=='min':
            value_cmap[value_cmap<0] = 0
        return value_cmap

#------------------------------------------------------------------------------
# Color stuff
#------------------------------------------------------------------------------
def hcl(c, degrees=False, gamma=3): # gamma is default
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
    h = h*180/np.pi if degrees else h/(2*np.pi) # normalize to 0-1 possibly
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

def cyclecolors(qcycle):
    """
    Just get the cycle colors.
    """
    # cycles = plt.get_cycles()
    cycles = mcolors.CYCLES
    cycle  = cycles.get(qcycle,None) # query
    if cycle is None:
        raise ValueError(f'Unknown cycle \"{qcycle}\". Options are: ' + ', '.join(cycles.keys()))
    return cycle

def colorshow(ncols=4, nbreak=12, minsat=0.1, cycle=False):
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
    for string in ['open','xkcd']:
        string = string or 'open'
        if string.lower()=='css':
            colors = mcolors.CSS4_COLORS # forget base colors, have them already; e.g. 'r' equals 'red'
        elif string.lower()=='xkcd':
            colors = mcolors.XKCD_SORTED
        elif string.lower()=='open':
            colors = mcolors.OPEN_COLORS
        else:
            raise ValueError(f"Unknown category \"{string}\".")
        # colors = {**mcolors.XKCD_SORTED}
        # Optionally cycle colors; they are not in the _colors_full_map list, so must add manually
        # But probably better to use cycleshow for this
        if cycle:
            seen = set() # trickery
            cycle_colors = rcParams['axes.prop_cycle'].by_key()['color']
            cycle_colors = [color for color in cycle_colors if not (color in seen or seen.add(color))] # trickery
            colors = {**colors, **{f'C{i}':v for i,v in enumerate(cycle_colors)}}
        # colors = {**mcolors.BASE_COLORS, **mcolors.XKCD_SORTED}
        # Sort colors by hue, saturation, value and name; then plot
        # Old, clunky method
        if string not in ('all','xkcd','open'):
            # Sort values
            by_hsv = sorted( (tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name) # an RGBtuple,name tuple
                for name,color in colors.items()) # sorts by first element in HSV tuple if different; then second; then third
            sorted_names = [name for hsv,name in by_hsv]
            # Create plot with mysterious methods
            figsize = (8, 5*(len(colors)/ncols)/40) # default is 5 inches tall when 40 colors present
            fig, ax = plt.subplots(figsize=figsize)
            nrows = len(sorted_names) // ncols + 1
            X, Y = fig.get_dpi() * fig.get_size_inches()
            h = Y/(nrows + 1)
            w = X/ncols
            for i,name in enumerate(sorted_names):
                col = i % ncols
                row = i // ncols
                y   = Y - (row * h) - h
                xi_line = w*(col + 0.05)
                xf_line = w*(col + 0.25)
                xi_text = w*(col + 0.3)
                displayname = name.split('xkcd:')[-1]
                n_names = len([n for n in sorted_names if n==displayname])
                if n_names>1:
                    print(f"Warning: \"{name}\" appears {n_names} times.")
                ax.text(xi_text, y, displayname, fontsize=h*0.8, ha='left', va='center')
                ax.hlines(y + h*0.1, xi_line, xf_line, color=colors[name], lw=h*0.6)
        # New method, that groups colors together by discrete range of hue then sorts by value
        else:
            # Fancy sort by grouping into hue categories, then sorting values; easy peasy
            # if string in ['open']: # group manually
            if string in ['open']:
                space = 0.5
                swatch = 1.5
                names = ["gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"]
                nrows, ncols = 10, len(names) # rows and columns
                sorted_names = [[name+str(i) for i in range(nrows)] for name in names]
            else: # just xkcd right now
                space = 1
                swatch = 1
                ncols = nbreak-1 # group by breakpoint
                colors_hsv = {k:tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(v))) for k,v in colors.items()}
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
        fig.savefig(f'{os.path.dirname(__file__)}/colors_{string}.pdf',
                bbox_inches='tight', format='pdf')
        figs += [fig]
    return figs

def cycleshow():
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
        seen = set()
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

