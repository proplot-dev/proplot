#!/usr/bin/env python3
from matplotlib import rcParams
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
from .rcmod import rc
from . import colortools as tools
from . import subplots # actually imports the function, since __init__ makes it global
_data = f'{os.path.dirname(__file__)}' # or parent, but that makes pip install distribution hard
_scales = {'rgb':(1,1,1), 'default':(360,100,100)}
_names  = {'rgb':('red', 'green', 'blue'),
            'hcl':('hue', 'chroma', 'luminance'),
            'hsl':('hue', 'saturation', 'luminance'),
            'hsv':('hue', 'saturation', 'value'),
            'hpl':('hue', 'partial sat', 'luminance')}

#------------------------------------------------------------------------------#
# Demo of channel values and colorspaces
#------------------------------------------------------------------------------#
def colorspace_breakdown(luminance=None, chroma=None, saturation=None, hue=None,
                         N=100, space='hcl'):
    """
    Generate hue-chroma, hue-luminance, and luminance-chroma
    cross-sections for the HCL, HSLuv, and HPLuv colorspaces.
    """
    # Dictionary
    hues = np.linspace(0, 360, 361)
    sats = np.linspace(0, 120, 120) # use 120 instead of 121, prevents annoying rough edge on HSL plot
    lums = np.linspace(0, 99.99, 101)
    chroma = saturation if saturation is not None else chroma
    if luminance is None and chroma is None and hue is None:
        luminance = 50
    if luminance is not None:
        hsl = np.concatenate((
            np.repeat(hues[:,None], len(sats), axis=1)[...,None],
            np.repeat(sats[None,:], len(hues), axis=0)[...,None],
            np.ones((len(hues), len(sats)))[...,None]*luminance,
            ), axis=2)
        suptitle = f'Hue-chroma cross-section for luminance {luminance}'
        xlabel, ylabel = 'hue', 'chroma'
        xloc, yloc = 60, 20
    elif chroma is not None:
        hsl = np.concatenate((
            np.repeat(hues[:,None], len(lums), axis=1)[...,None],
            np.ones((len(hues), len(lums)))[...,None]*chroma,
            np.repeat(lums[None,:], len(hues), axis=0)[...,None],
            ), axis=2)
        suptitle = f'Hue-luminance cross-section for chroma {chroma}'
        xlabel, ylabel = 'hue', 'luminance'
        xloc, yloc = 60, 20
    elif hue is not None:
        hsl = np.concatenate((
            np.ones((len(lums), len(sats)))[...,None]*hue,
            np.repeat(sats[None,:], len(lums), axis=0)[...,None],
            np.repeat(lums[:,None], len(sats), axis=1)[...,None],
            ), axis=2)
        suptitle = 'Luminance-chroma cross-section'
        xlabel, ylabel = 'luminance', 'chroma'
        xloc, yloc = 20, 20

    # Make figure, with hatching indiatinc invalid values
    # Note we invert the x-y ordering for imshow
    # rc['facehatch'] = '....'
    rc['facecolor'] = 'k'
    # rc['facehatch'] = 'xxx'
    f, axs = subplots(ncols=3, bottomlegends=True, rightcolorbar=True,
                        span=0, share=0, wspace=0.6, axwidth=2.5,
                        bottom=0, left=0, right=0,
                        aspect=1, tight=True)
    for i,(ax,space) in enumerate(zip(axs,('hcl','hsl','hpl'))):
        rgba = np.ones((*hsl.shape[:2][::-1], 4)) # RGBA
        for j in range(hsl.shape[0]):
            for k in range(hsl.shape[1]):
                rgb_jk = tools.to_rgb(hsl[j,k,:].flat, space)
                # rgba[k,j,:3] = np.clip(rgb_jk, 0, 1)
                if not all(0 <= c <= 1 for c in rgb_jk):
                    rgba[k,j,3] = 0 # transparent cell
                else:
                    rgba[k,j,:3] = rgb_jk
        ax.imshow(rgba, origin='lower', aspect='auto')
        ax.format(xlabel=xlabel, ylabel=ylabel, suptitle=suptitle,
                  grid=False, tickminor=False,
                  xlocator=xloc, ylocator=yloc,
                  title=space.upper(), title_kw={'weight':'bold'})
    return f

def cmap_breakdown(cmap, N=100, space='hcl'):
    """
    Handy visualization of how a `~proplot.colortools.PerceptuallyUniformColormap`
    varies in the HCL, HSLuv, and HPLuv colorspaces.
    """
    # Figure
    f, axs = subplots(ncols=4, bottomlegends=True, rightcolorbar=True,
                           span=0, sharey=1, wspace=0.5,
                           bottom=0.4, axwidth=2, aspect=1, tight=True)
    x = np.linspace(0, 1, N)
    cmap = tools.colormap(cmap, N=N) # arbitrary cmap argument
    cmap._init()
    name = cmap.name
    for j,(ax,space) in enumerate(zip(axs,('hcl','hsl','hpl','rgb'))):
        # Get RGB table, unclipped
        hs = []
        if hasattr(cmap, 'space'):
            lut = cmap._lut_hsl[:,:3].copy()
            for i in range(len(lut)):
                lut[i,:] = tools.to_rgb(lut[i,:], cmap.space)
        else:
            lut = cmap._lut[:,:3].copy()
        # Convert RGB to space
        for i in range(len(lut)):
            lut[i,:] = tools.to_xyz(lut[i,:], space=space)
        scale  = _scales.get(space, _scales['default'])
        labels = _names[space]
        # Draw line, add legend
        colors = ['C1', 'C2', 'C0'] # corresponds with RGB roughly
        m = 0
        for i,label in enumerate(labels):
            y = lut[:-2,i]/scale[i]
            y = np.clip(y, 0, 5)
            h, = ax.plot(x, y, color=colors[i], lw=2, label=label)
            m = max(m, max(y))
            hs += [h]
        f.bottompanel[j].legend(hs)
        ax.axhline(1, color='gray7', dashes=(1.5, 2.5), alpha=0.8, zorder=0, lw=2)
        ax.format(title=space.upper(), titlepos='oc', ylim=(0-0.1, m + 0.1))
    # Draw colorbar
    with np.errstate(all='ignore'):
        m = ax.contourf([[np.nan,np.nan],[np.nan,np.nan]], levels=100, cmap=cmap)
    f.rightpanel.colorbar(m, clocator='none', cformatter='none', clabel=f'{name} colors')
    locator = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 10]
    axs.format(suptitle=f'{name} colormap breakdown', ylim=None, ytickminor=False,
              yscale=('cutoff', 4, 1), ylocator=locator, # progress 10x faster above x=1
              xlabel='position', ylabel='scaled channel value')

#------------------------------------------------------------------------------#
# Reference tables for colors, colormaps, cycles
#------------------------------------------------------------------------------#
def color_show(groups=None, ncols=4, nbreak=12, minsat=0.2):
    """Visualize all possible named colors. Adapted from `this example
    <https://matplotlib.org/examples/color/named_colors.html>`_."""
    # Get colors explicitly defined in _colors_full_map, or the default
    # components of that map (see soure code; is just a dictionary wrapper
    # on some simple lists)
    figs = []
    scale = (360, 100, 100)
    groups = groups or [['xkcd','crayons']]
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
                color_dict.update(tools.colorlist[name]) # add category dictionary

        # Group colors together by discrete range of hue, then sort by value
        # For opencolors this is not necessary
        if 'open' in group:
            # Sorted color columns and plot settings
            wscale = 0.5
            swatch = 1.5
            names = ['red', 'pink', 'grape', 'violet', 'indigo', 'blue', 'cyan', 'teal', 'green', 'lime', 'yellow', 'orange', 'gray']
            nrows, ncols = 10, len(names) # rows and columns
            plot_names = [[name + str(i) for i in range(nrows)] for name in names]
            nrows = nrows*2
            ncols = (ncols+1)//2
            plot_names = np.array(plot_names, order='C')
            plot_names.resize((ncols, nrows))
            plot_names = plot_names.tolist()
        else:
            # For other palettes this is necessary
            # Get colors in perceptally uniform space
            # Then will group based on hue thresholds
            wscale = 1
            swatch = 1
            colors_hsl = {
                key: [c/s for c,s in zip(tools.to_xyz(value, tools._distinct_colors_space), scale)]
                for key,value in color_dict.items()
                }

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
            # ncols = nbreak-1 # allow custom setting
            names = [i for sublist in plot_names for i in sublist]
            plot_names = [[]]
            nrows = len(names)//ncols+1
            for i,name in enumerate(names):
                if ((i + 1) % nrows)==0:
                    plot_names.append([]) # add new empty list
                plot_names[-1].append(name)

        # Create plot by iterating over columns 
        # Easy peasy. And put 40 colors in a column
        fig, ax = subplots(width=8*wscale*(ncols/4),
                           height=5*(nrows/40),
                           left=0, right=0, top=0, bottom=0,
                           tight=False)
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
                ax.hlines(y_line, xi_line, xf_line, color=color_dict[name], lw=hsep*0.6)
        # Save
        # fig.save(f'{_data}/colors/colors_{"-".join(group)}.pdf', format='pdf', transparent=False)
        ax.format(xlim=(0,X), ylim=(0,Y))
        ax.set_axis_off()
        figs += [fig]
    return figs

def cycle_show():
    """Show off the different color cycles."""
    # Get the list of cycles
    _cycles = {name:mcm.cmap_d[name].colors for name in tools.cycles}
    _cycles = {name:_cycles[name] for name in sorted(_cycles.keys())}
    nrows = len(_cycles)//2 + len(_cycles)%2
    # Create plot
    state = np.random.RandomState(528)
    fig, axs = subplots(width=6, wspace=0.05, hspace=0.25,
                        sharey=False, sharex=False,
                        aspect=2, ncols=2, nrows=nrows)
    for i,(ax,(key,cycle)) in enumerate(zip(axs, _cycles.items())):
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
    if len(_cycles)%2==1:
        axs[-1].set_visible(False)
    # Save
    # fig.save(f'{_data}/colors/cycles.pdf', format='pdf')
    return fig

# def cmap_show(N=31):
def cmap_show(N=129):
    """Plots all current colormaps, along with their catgories. Adapted from
    `this example <http://matplotlib.org/examples/color/colormaps_reference.html>`_."""
    # Have colormaps separated into categories
    cmaps_reg = [name for name in mcm.cmap_d.keys() if name not in ('vega', 'greys')
            and isinstance(mcm.cmap_d[name], mcolors.LinearSegmentedColormap)]

    # Detect unknown/manually created colormaps, and filter out
    # colormaps belonging to certain section
    categories = {cat:names for cat,names in tools._cmap_categories.items()
                        if cat not in tools._cmap_categories_delete}
    categories_reg = {cat:[name for name in names if name.lower() in cmaps_reg]
                        for cat,names in categories.items()}
    cmaps_delete  = [name.lower() for cat,names in tools._cmap_categories.items() for name in names
                        if cat in tools._cmap_categories_delete]
    cmaps_reg_known   = [name.lower() for cat,names in categories.items() for name in names
                        if name.lower() in cmaps_reg]
    cmaps_missing = [name.lower() for cat,names in categories.items() for name in names
                        if name.lower() not in cmaps_reg]
    cmaps_unknown  = [name for name in cmaps_reg
                        if name not in cmaps_reg_known and name not in cmaps_delete]
    if cmaps_missing:
        print(f'Missing colormaps: {", ".join(cmaps_missing)}')
    if cmaps_delete:
        print(f'Deleted colormaps: {", ".join(cmaps_delete)}')
    if cmaps_unknown:
        print(f'New colormaps: {", ".join(cmaps_unknown)}')

    # Array for producing visualization with imshow
    a = np.linspace(0, 1, 257).reshape(1,-1)
    a = np.vstack((a,a))
    # Figure
    extra = 1 # number of axes-widths to allocate for titles
    nmaps = len(cmaps_reg_known) + len(cmaps_unknown) + len(categories)*extra
    fig, axs = subplots(nrows=nmaps, axwidth=4.5, axheight=0.23, span=False, share=False, hspace=0.07)
    iax = -1
    ntitles, nplots = 0, 0 # for deciding which axes to plot in
    for cat,names in categories_reg.items():
        # Space for title
        ntitles += extra # two axes-widths
        for imap,name in enumerate(names):
            # Checks
            iax += 1
            if imap + ntitles + nplots > nmaps:
                ax.invisible()
                break
            ax = axs[iax]
            if imap==0:
                iax += 1
                ax.invisible()
                ax = axs[iax]
            if name not in mcm.cmap_d or name.lower() not in cmaps_reg: # i.e. the expected builtin colormap is missing
                ax.invisible() # empty space
                continue
            # Draw map
            ax.imshow(a, cmap=name, origin='lower', aspect='auto', levels=N)
            ax.format(ylabel=name, ylabel_kw={'rotation':0, 'ha':'right', 'va':'center'},
                      xticks='none',  yticks='none', # no ticks
                      xloc='neither', yloc='neither', # no spines
                      title=(cat if imap==0 else None))
        # Space for plots
        nplots += len(names)
    # Save
    # fig.save(f'{_data}/cmaps/colormaps.pdf')
    return fig

