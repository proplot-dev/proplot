#------------------------------------------------------------------------------
# All
#------------------------------------------------------------------------------
__all__ = [
    'subplots', 'format', 'basemap', # basic setup
    'Normalize', # custom colormap normalization
    ]

#------------------------------------------------------------------------------
# Imports, all
#------------------------------------------------------------------------------
from . import utils
import os # need for saving figures
import matplotlib as mpl
import mpl_toolkits.basemap as basemap
import matplotlib.pyplot as plt # but should already be imported
import numpy as np
import pandas as pd
import cartopy.crs as ccrs # crs stands for "coordinate reference system", leading c is "cartopy"
import cartopy.mpl as cmpl

#------------------------------------------------------------------------------
# Definitions
#------------------------------------------------------------------------------
def format(ax, filename=None, desktop=True, # file name, and whether to put in desktop folder
        size=None, fix=False, # figure size, and whether to fix aspect ratio
        legend=True, objects=None, labels=None, legend_ax=None, legend_kw=dict(), # legend options; if declared and have no objects with 'label', does nothing
        colorbar=False, mappable=None, colorbar_ax=None, colorbar_kw=dict(), # colorbar options
        xgrid=True, ygrid=True, xgridminor=True, ygridminor=True, # grid options (add to this)
        xtickminor=False, ytickminor=False, # turning on/off minor ticks
        xtickpos='bottom', ytickpos='left', xopposite=True, yopposite=True, # location of ticks, which spines to draw
        xlim=None, ylim=None, xscale='linear', yscale='linear', xreverse=False, yreverse=False, # special properties
        xdates=False, ydates=False, # whether to format axis labels as long datetime strings
        title=None, suptitle=None, abc=False, # label subplots with a), b), c)
        xlabel=None, ylabel=None, clabel=None, # axis/colorbar labels
        xlocator=None, xminorlocator=None, clocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
        xformatter=None, cformatter=None, yformatter=None, # formatter
        color='k', linewidth=0.7, # for **all** line objects (legend/axis/colorbar)
        ctick_kw=dict(), xtick_kw=dict(), ytick_kw=dict(), # all should be special tick properties, using set_tick_params; also takes axis_kw color
        ticklabels_kw=dict(), label_kw=dict(), title_kw=dict(), suptitle_kw=dict(), abc_kw=dict(), # all are passed to FontProperties class; also takes axis_kw color
        grid_kw=dict(), gridminor_kw=dict() # for axis properties, all are line/patch-type properties
        ):
    """
    TODO: Add options for labels spanning multiple subplots.
    Possibly forget some keyword groups; maybe just take 'color' and 'linewidth' for 
    axis limits and fonts, with grid.
    
    TODO: Write module for custom lat/lon annotations... should be able to be done
    with .annotate specifying lon/lat coordinates, I believe.
    
    NOTE: possible date axes objects are TimeStamp (pandas), np.datetime64, DateTimeIndex; 
    can fix with fig.autofmt_xdate() or manually set options; uses ax.is_last_row()
    or ax.is_first_column(), so should look that up...problem is there is no 
    autofmt_ydate(), so really should implement my own version of this...

    NOTE: sharex and sharey only affect *ticklocations*/*ticklabels*/*tickformatters*/*limits*
    by making axes correspond to each other. Does *not* affect axis labels, unless you used
    my.subplots (which set them invisible); this function also makes them invisible.
    """
    # Allow applying identical settings to *multiple* axes at once
    try: iter(ax)
    except TypeError: ax = np.array([ax])
    else: ax = np.array(ax) # if already array, this does no harm; like list([1,2,3])

    # Override defaults for legend/colorbar, labels, ticklabels, title, grid/gridminor
    # ...can be modified
    legend_kw = utils.dict(dict(framealpha=1, fancybox=False), **legend_kw)
    ticklabels_kw = utils.dict(dict(size=8, weight='normal', color=color), **ticklabels_kw)
    label_kw = utils.dict(dict(size=8, weight='normal', color=color), **label_kw)
    title_kw = utils.dict(dict(size=8, weight='bold', color=color), **title_kw)
    abc_kw = utils.dict(dict(size=8, weight='bold', color=color), **abc_kw)
    colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=False), **colorbar_kw)
    grid_kw = utils.dict(dict(linestyle='-', linewidth=0.8, color=color, alpha=0.1), **grid_kw)
    gridminor_kw = utils.dict(dict(linestyle=':', linewidth=0.5, color=color, alpha=0.05), **gridminor_kw)
    xtick_kw = utils.dict(dict(length=5, direction='in', width=linewidth, color=color), **xtick_kw)
    ytick_kw = utils.dict(dict(length=3, direction='in', width=linewidth, color=color), **ytick_kw)
    ctick_kw = utils.dict(dict(length=5, direction='out', width=linewidth, color=color), **ctick_kw)
    # ...are not modified (depend on others)
    xtickminor_kw = utils.dict(xtick_kw, length=xtick_kw.length/2)
    ytickminor_kw = utils.dict(ytick_kw, length=ytick_kw.length/2)
    coastlines_kw = utils.dict(linewidth=linewidth, color=color)
    gridlines_kw = utils.dict(linestyle=':', linewidth=linewidth, color=color, alpha=0.2,
            draw_labels=False, crs=ccrs.PlateCarree())
    # ...misc
    abc_convert = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e', 5:'f', 6:'g', 7:'h', 8:'j', 9:'i', 10:'k', 11:'l', 12:'m'}

    # Colorbar/legend checks
    if colorbar and mappable is None: # steal space from existing axis
        raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
    if legend and objects is None and legend_ax is not None:
        raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')

    # Set properties
    for i,a in enumerate(ax.flat): # .flat gives row-major indexing, just like matlab for plots
        skip = False # to be used later
        # Legend axes setup
        # ...fist, make object
        leg = None
        if legend and legend_ax is None:
            leg = a.legend(**legend_kw)
        elif legend and (legend_ax % ax.size)==i: # e.g. can input -1, or axes.Axes object
            a.set_frame_on(False)
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)
            leg = a.legend(handles=objects, loc='center', **legend_kw)
            skip = True
        # ...next, set properties
        if leg is not None:
            # ...check for certain properties only
            leg.get_frame().set_edgecolor(color)
            leg.get_frame().set_linewidth(linewidth)
            # ...and apply text properties, identical to ticklabel properties
            plt.setp(leg.get_texts(), **ticklabels_kw) # color should be in there as font property

        # Colorbar axes setup
        # ...first, make object input
        special_kw = None
        if colorbar and colorbar_ax is None:
            special_kw = dict(ax=a) # steals space form axis 'a', uses make_axes() to build new one
        elif colorbar and (colorbar_ax % ax.size)==i:
            skip = True
            special_kw = dict(cax=a, use_gridspec=True) # uses space affored by this entire axis
        # ...next, set properties
        if special_kw is not None:
            # ...create object
            cb = a.get_figure().colorbar(mappable, **special_kw, **colorbar_kw)
            # ...first, make edges/dividers consistent with axis edges
            if cb.dividers is not None: plt.setp(cb.dividers, color=color, linewidth=linewidth) # change dividers
            cb.outline.set_edgecolor(color)
            cb.outline.set_linewidth(linewidth)
            # ...next, ticks formatters/locators (different from axis_setup version)
            if isinstance(clocator, mpl.ticker.Locator):
                cb.locator = clocator
            elif clocator is not None:
                try: iter(clocator)
                except TypeError: clocator = utils.arange(cb.vmin, cb.vmax, clocator)
                cb.locator = mpl.ticker.FixedLocator(clocator)
            if isinstance(cformatter, mpl.ticker.Formatter):
                cb.formatter = cformatter
            elif cformatter is not None:
                if isinstance(cformatter, str):
                    cb.formatter = mpl.ticker.FormatStrFormatter(cformatter) # %-style, numbers
                else:
                    cb.formatter = mpl.ticker.FixedFormatter(cformatter)
            # ...and update (use this instead of update_bruteforce)
            cb.update_ticks() # updates formatters/locators; colorbar uses weird combo of its own patch objects, and underlying axis
            # ...finally, the ticks/ticklabels/labels basic properties
            if cb.orientation=='horizontal':
                axis = cb.ax.xaxis
            else:
                axis = cb.ax.yaxis
            plt.setp(axis.get_ticklabels(), **ticklabels_kw)
            axis.set_tick_params(which='both', **ctick_kw)
            if clabel is not None: axis.set_label_text(clabel, **label_kw)
            # ...how can i do this
            # ...problem is axis has vmin, vmax
            # ...and next, tick stuff (hard to deal with lines separating colors, vs. tick marks)
        
        # Projection axes setup
        if hasattr(a, 'projection'): # is then a subclass generated from Cartopy
            skip = True
            a.coastlines(**coastlines_kw)
            gl = a.gridlines(**gridlines_kw)
            if gridlines_kw.draw_labels:
                gl.xlabels_top = False
                gl.ylabels_left = False
                gl.xlocator = mpl.ticker.FixedLocator(utils.arange(-180, 180, 60))
                gl.ylocator = mpl.ticker.FixedLocator(utils.arange(-90, 90, 30))
                gl.xformatter = cmpl.gridliner.LONGITUDE_FORMATTER
                gl.yformatter = cmpl.gridliner.LATITUDE_FORMATTER
                gl.xlabel_style = gl.ylabel_style = ticklabels_kw
                # f.canvas.draw()
        
        # Skip, if axes was either meant to be filled with legend or colorbar, or was projection
        if skip:
            continue

        # Control visibility, if axis has sister shared axis
        if a._sharex is not None:
            plt.setp(a.xaxis.get_ticklabels(), visible=False) 
            a.xaxis.get_label().set_visible(False)
        if a._sharey is not None:
            plt.setp(a.yaxis.get_ticklabels(), visible=False) 
            a.yaxis.get_label().set_visible(False)

        # Axis properties, that must be set in unison or have special properties
        # ...spine line properties, tick line properties (should be identical between axes)
        plt.setp(tuple(a.spines.values()), color=color, linewidth=linewidth) # axis properties; universal
        # ...tick alliegance
        for axis, tickpos in zip((a.xaxis, a.yaxis), (xtickpos, ytickpos)):
            if tickpos is not None:
                axis.set_ticks_position(tickpos) # both/left/right
        # ...spines visibility
        for axis, sides, opposite in zip((a.xaxis, a.yaxis), (('bottom','top'), ('left','right')), (xopposite,yopposite)):
            for side in sides:
                if axis.get_ticks_position() in ('both',side) or opposite:
                    a.spines[side].set_visible(True) # view top if ticks there, or turned it on
                elif not opposite:
                    a.spines[side].set_visible(False)
        # ...axis scaling
        a.set_xscale(xscale)
        a.set_yscale(yscale)
        # ...axis limits
        if xlim is None: xlim = a.get_xlim()
        if ylim is None: ylim = a.get_ylim()
        if xreverse: xlim = xlim[::-1]
        if yreverse: ylim = ylim[::-1]
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        # ...and grid stuff
        if xgrid:
            a.xaxis.grid(which='major', **grid_kw)
        if ygrid:
            a.yaxis.grid(which='major', **grid_kw)
        if xgridminor and xtickminor:
            a.xaxis.grid(which='minor', **gridminor_kw)
        if ygridminor and ytickminor:
            a.yaxis.grid(which='minor', **gridminor_kw)

        # X/Y Axis properties
        for axis, label, dates, gridminor, tickminor, tickminorlocator, tickminor_kw, grid, ticklocator, tickformatter, tick_kw in zip(
                (a.xaxis, a.yaxis), (xlabel, ylabel), (xdates, ydates), # other stuff
                (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), (xtickminor_kw, ytickminor_kw), # minor ticks
                (xgrid, ygrid), (xlocator, ylocator), (xformatter, yformatter), (xtick_kw, ytick_kw) # major ticks
                ):
            # ...first, major tick locators (should not affect limits)
            lim = axis.get_view_interval() # to be used, automatically
            if isinstance(ticklocator, mpl.ticker.Locator):
                axis.set_major_locator(ticklocator)
            elif ticklocator is not None:
                try: iter(ticklocator)
                except TypeError: ticklocator = utils.arange(lim[0], lim[1], ticklocator)
                # axis.set_ticks(ticklocator)
                axis.set_major_locator(mpl.ticker.FixedLocator(ticklocator))
            # ...next, minor tick locator
            if not tickminor:
                axis.set_minor_locator(mpl.ticker.NullLocator())
            elif isinstance(tickminorlocator, mpl.ticker.Locator):
                axis.set_minor_locator(tickminorlocator)
            elif tickminorlocator is not None:
                try: iter(tickminorlocator)
                except TypeError: tickminorlocator = utils.arange(lim[0], lim[1], tickminorlocator)
                axis.set_minor_locator(mpl.ticker.FixedLocator(tickminorlocator))
            # ...next, tick formatters (enforce Null, always, for minor ticks)
            axis.set_minor_formatter(mpl.ticker.NullFormatter())
            if isinstance(tickformatter, mpl.ticker.Formatter): # use isinstance, because r"stuff" or f"stuff" **do not work**
                axis.set_major_formatter(tickformatter)
            elif tickformatter is not None:
                if isinstance(tickformatter, str):
                    if dates: axis.set_major_formatter(mpl.dates.DateFormatter(tickformatter)) # %-style, dates
                    else: axis.set_major_formatter(mpl.ticker.FormatStrFormatter(tickformatter)) # %-style, numbers
                else:
                    # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
                    axis.set_major_formatter(mpl.ticker.FixedFormatter(tickformatter)) # list of strings
            plt.setp(axis.get_ticklabels(), **ticklabels_kw)
            # ...and finally, label
            if label is not None: axis.set_label_text(label)
            plt.setp(axis.get_label(), **label_kw)
            # ...and tick properties
            axis.set_tick_params(which='major', **tick_kw)
            axis.set_tick_params(which='minor', **tickminor_kw) # have length
        
        # Title, and abc'ing, and suptitle
        if title is not None:
            a.set_title(title[i], **title_kw)
        if abc:
            a.annotate(abc_convert[i] + ')', xy=(0,0), xytext=(0,1.03),
                    xycoords='axes fraction', textcoords='axes fraction',
                    verticalalignment='baseline', horizontalalignment='left',
                    **title_kw)
    
    # Options for entire figure
    if suptitle is not None:
        a.get_figure().suptitle(suptitle, **title_kw)
    if size is not None:
        if type(size) is str:
            pass # allow for presets for specific journals
        a.get_figure().set_size_inches(size)
    if fix:
        a.get_figure().set_size_inches() # determine appropriate size, given 
            # ...number of rows/columns and subplot aspect ratios; should work
            # ...because even colorbars/legends have their own axes

    # And save, optionally (into figures folder)
    if filename is not None:
        if '.' not in filename:
            filename = filename + '.pdf'
        if desktop: 
            filename = os.environ['HOME'] + '/Desktop/Figures/' + filename
        a.get_figure().savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.05, dpi=300,
                transparent=True, frameon=False)
            # pad_inches is padding around bbox, and frameon makes figure window transparent
            # even though saving as pdf, think dpi will affect any embedded raster objects
            
def subplots(array=None, nrows=1, ncols=1, sharex=False, sharey=False, # for sharing x, y axes
        fixedaspect=False, aspect=None, height=None, width=7, # for controlling aspect ratio; default is control for width
        panelx=False, panely=False, panelx_ratio=0.05, panely_ratio=0.05, # add panels for legends, colorbars, etc.
        hspace=0.2, wspace=0.2, height_ratios=None, width_ratios=None, # options for gridspec (instead of gridpsec_kw)
        figsize=None, # options for plt.figure (right now, just figure size)
        projection=None, **projection_kw): # options for add_subplot (projections, etc.; 
            # and accept general projection kwargs, because they very a lot between choices)
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want a smallish panel
    to plot data, just use width/height ratios.
    """
    # First things first, interpret rarray input and ratios
    # ...array stuff
    if array is None:
        array = np.reshape(np.arange(1,nrows*ncols+1), (nrows, ncols))
    array = np.array(array)
    nrows, ncols = array.shape[0], array.shape[1]
    # ...width/height ratios
    if width_ratios is None:
        width_ratios = np.ones(ncols)/ncols
    else:
        width_ratios = np.array(width_ratios)/sum(width_ratios)
    if height_ratios is None:
        height_ratios = np.ones(nrows)/nrows
    else:
        height_ratios = np.array(height_ratios)/sum(height_ratios)
    # ...figure size default
    if figsize is None:
        figsize = (width, width*nrows/ncols)
    
    # Create dictionary converting projections to classes
    if projection is not None:
        # ...make axis
        crs_dict = dict(
                rectilinear=ccrs.PlateCarree,
                pcarree=ccrs.PlateCarree,
                platecarree=ccrs.PlateCarree, 
                robinson=ccrs.Robinson,
                stereo=ccrs.Stereographic,
                stereographic=ccrs.Stereographic,
                moll=ccrs.Mollweide,
                mollweide=ccrs.Mollweide,
                )
        projection = crs_dict[projection](**projection_kw)
        # ...and aspect ratio considerations
        fixedaspect = True
        x, y = projection.x_limits, projection.y_limits # in projection coordinates
        aspect = (y[1]-y[0])/(x[1]-x[0])
        aspect = 1

    # Panels for colorbar/legend; adjust input array, and width/height ratios
    if panely: # right panel
        width_ratios = np.append(width_ratios, panelx_ratio)/(1+panelx_ratio)
        array = np.concatenate((array, np.ones((nrows, 1))*(array.max()+1)), axis=1)
        ncols += 1
    if panelx: # bottom panel
        height_ratios = np.append(height_ratios, panely_ratio)/(1+panely_ratio)
        array = np.concatenate((array, np.ones((1, ncols))*(array.max()+1)), axis=0)
        nrows += 1
    
    # Fixed aspect, infer required figure size (works because hspace/wspace are in terms of the average SubplotSpec 
    # [indexed GridSpec object], not of the realized axes drawn from gridspec objects [try wspace=0.2 for gs with 
    # cols=4, axes gs[0,:2] and gs[0,2:] vs. gs with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ])
    if fixedaspect and height is not None:
        # ...set figure aspect ratio, fixing height
        width = height*(nrows + (nrows-1)*wspace) / (aspect*(ncols + (ncols-1)*hspace))
        figsize = (width, height)
    elif fixedaspect:
        # ...set figure aspect ratio, fixing width
        height = width*(aspect*(ncols + (ncols-1)*hspace)) / (nrows + (nrows-1)*wspace)
        figsize = (width, height)
    # ...also, don't worry about left/bottom/top/right, because we will save figure
    # with tight bounding box anyway... or not

    # Set up figure and gridspec object
    fig = plt.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(nrows=array.shape[0], ncols=array.shape[1], 
            wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)

    # Some iniital/important stuff for getting shared axes
    if panely: # gets priority
        array = array[:,:-1]
    if panelx:
        array = array[:-1,:] # for testing purposes
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # ...0 or -1 stands for 'empty'; could also use NaN
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    num_axes = len(axes_ids)
    
    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # ...get shared x
    if sharex:
        xgroups_share, xgroups, grouped = [], [], []
        for i in range(num_axes):
            matches = np.all(xrange[i:i+1,:]==xrange, axis=1) # broadcasting rules work here
            group = np.where(matches)[0]
            if i not in grouped and group.size>1:
                xgroups += [group] # add ndarray of ids to list
                xgroups_share += [group[np.argmax(yrange[group, 1])]] # bottom-most axis with shared x, for group
            grouped += [*group] # bookkeeping; record ids that have been grouped already
    # ...get shared y
    if sharey:
        ygroups_share, ygroups, grouped = [], [], []
        for i in range(num_axes):
            matches = np.all(yrange[i:i+1,:]==yrange, axis=1) # broadcasting rules work here
            group = np.where(matches)[0]
            if i not in grouped and group.size>1:
                ygroups += [group] # add ndarray of ids to list
                ygroups_share += [group[np.argmin(xrange[group, 0])]] # left-most axis with shared y, for group
            grouped += [*group] # bookkeeping; record ids that have been grouped already
    
    # Make shared axes
    ax = np.array(num_axes*[None]) # is array with dtype 'object'
    if sharex:
        for i in xgroups_share:
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], projection=projection)
    if sharey:
        for i in ygroups_share:
            if ax[i] is None:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], projection=projection)

    # Make dependent axes
    # ...have to manually insert the correct axis object for sharex, sharey argument on creation with add_subplot
    # print(num_axes, array, '\nxrange:', xrange, '\nyrange:', yrange, '\nxgroups, and shared:', xgroups, '\n', xgroups_share, '\nygroups, and shared:', ygroups, '\n', ygroups_share)
    for i in range(num_axes):
        if ax[i] is None: # i.e. it's not one of the axes on which others are dependent
            # ...detect if we want to share this axis with another
            sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
            if sharex:
                igroup = np.where([i in g for g in xgroups])[0] # np.where works on lists
                if igroup.size==1:
                    sharex_ax = ax[xgroups_share[igroup[0]]]
                    if sharex_ax is None:
                        raise ValueError('Something went wrong; shared x axes was not already drawn.')
                elif igroup.size>1:
                    raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)
            if sharey:
                igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
                if igroup.size==1:
                    sharey_ax = ax[ygroups_share[igroup[0]]]
                    if sharey_ax is None:
                        raise ValueError('Something went wrong; shared x axes was not already drawn.')
                elif igroup.size>1:
                    raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)
            # ...draw axis, and add to list
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], sharex=sharex_ax, sharey=sharey_ax,
                    projection=projection)
            # ...hide tick labels (not default behavior for manual sharex, sharey use)
            if sharex_ax is not None:
                plt.setp(ax[i].xaxis.get_ticklabels(), visible=False) 
                ax[i].xaxis.get_label().set_visible(False)
                # setp acts for iterable set of matplotlib objects; this is not black magic
                # these hold, because setting ticklabels or labels with higher-level commands
                # just changes the **text**, not underlying properties... pretty neat!
            if sharey_ax is not None:
                plt.setp(ax[i].yaxis.get_ticklabels(), visible=False)
                ax[i].yaxis.get_label().set_visible(False)
                
    # Extra stuff
    # ...get transform objcts for each axes
    if projection is not None:
        # ...really convenient; axes have their own projection, but can plot data gathered from graticule/mesh 
        # with **different** projection on top; in our case, data is almost always lon/lat rectinlinear
        transform = np.array(num_axes*[None]) # is array with dtype 'object'
        for i in range(num_axes):
            ax[i].transform = ccrs.PlateCarree()._as_mpl_transform(ax[i])
            # ...have to use _as_mpl_transform, or annotations won't work
    # ...set property cycles, for line plots, and other axes properties
    for i in range(num_axes):
        ax[i].set_prop_cycle(mpl.cycler(
            color=(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']*4),
            linestyle=(['-']*10 + ['--']*10 + [':']*10 + ['-.']*10)
            ))
    # ...create panel axes (at end)
    if panely:
        a = fig.add_subplot(gs[:,-1])
        ax = np.append(ax, a)
    if panelx:
        a = fig.add_subplot(gs[-1,:])
        ax = np.append(ax, a)

    # Finally, return result
    return fig, ax

def cmapload(name):
    """
    Load colormaps from RGB file using pandas.
    Call this function whenever startying up matplotlib.
    """
    data = pd.read_table(name)
    has_neutral, has_cutoff = False, False
    if type(name) is str:
        if name=='bw':
            cmap = np.tile(utils.arange(0, 1, 0.1)[np.newaxis,:], (1,3))
        # Rainbows
        elif name=='rainbows':
            cmap = pd.read_table(cmaploc + 'WhiteBlueGreenYellowRed.rgb', sep='\s+', header=None)
            cmap = cmap/255
        # Blue-reds: hard
        elif name=='bwo': # hot reddish-orange to blue, white in middle
            cmap = pd.read_table(cmaploc + 'BlueWhiteOrangeRed.rgb', sep='\s+', header=None)
            cmap = cmap/255 # has white in middle, blue to light blue to white, yellow, red
            has_neutral = True
        elif name=='br':
            cmap = pd.read_table(cmaploc + 'BlRe.rgb', sep='\s+', header=None)
            cmap = cmap/255 # red-blue very saturated, with sharp transition
            has_cutoff = True
        elif name=='brsoft':
            cmap = pd.read_table(cmaploc + 'BlueRed.rgb', skiprows=1, sep='\s+', header=None)
            cmap = cmap/255 # NEW
            has_cutoff = True
        elif name=='bo':
            cmap = pd.read_table(cmaploc + 'BlueYellowRed.rgb', sep='\s+', header=None)
            cmap = cmap/255 # blue to orange/red, with SHARP transition in middle
            has_cutoff = True
        elif name=='wbrw':
            cmap = pd.read_table(cmaploc + 'WhBlReWh.rgb', skiprows=1, sep='\s+', header=None)
            cmap = cmap/255 # NEW
            has_cutoff = True
        # Land
        elif name=='terrain':
           cmap = pd.read_table(cmaploc + 'OceanLakeLandSnow.rgb', sep='\s+', header=None)
           cmap = cmap[2:,:]/255 # ignore the ocean/lake colors
        elif name=='night':
            cmap = pd.read_table(cmaploc + 'GMT_nighttime.rgb', skiprows=2, sep='\s+', header=None) # NEW
            has_cutoff = True
        # Dry-wets
        elif name=='drywet':
            cmap = pd.read_table(cmaploc + 'GMT_drywet.rgb', skiprows=2, sep='\s+', header=None)
        elif name=='wgb': # very few levels
            cmap = pd.read_table(cmaploc + 'CBR_wet.rgb', skiprows=2, sep='\s+', header=None)
            cmap = cmap/255
        elif name=='wgbv': # very few levels
            cmap = pd.read_table(cmaploc + 'precip_11lev.rgb', skiprows=2, sep='\s+', header=None)
            cmap = cmap/255
        elif name=='bwgb': # brown white green blue
            cmap = pd.read_table(cmaploc + 'precip_diff_12lev.rgb', skiprows=2, sep='\s+', header=None)
            cmap = cmap/255
            has_neutral = True
        # Aliases for existing maps
        # ...put code here
        # ...put code here
        else:
            raise ValueError('Unknown colormap string.')
    elif type(name) is np.ndarray:
        # User input their own colormap, and wants to interpolate/modify it
        cmap = name
    else:
        raise ValueError('Unknown input type; must be string or numpy array.')
    if cmap.shape[1] not in (3, 4):
        raise ValueError('Bad colormap input; must be RGB or RGBA.')
    # And create object
    cmap_dict = utils.dict(red=cmap[:,0], green=cmap[:,1], blue=cmap[:,2])
    if cmap.shape[1]==4: cmap_dict.alpha = cmap[:,3]
    cmap_object = mpl.colors.LinearSegmentedColormap(name, cmap_dict, cmap.shape[0])
    plt.register_cmap(cmap=cmap_object)
    return cmap_object

# def basemap(ax, projection='moll', lon_0=0):
#     """
#     Automatic construction of basemap for each axes.
#     """
#     # Convert to array
#     try:
#         iter(ax)
#     except TypeError:
#         ax = np.array([ax])
#     else:
#         ax = np.array(ax) # if already array, this does no harm; like list([1,2,3])
#     # Iterate
#     for a in ax:
#         # draw map, and some other useful stuff
#         m = basemap.basemap.Basemap(ax=a,projection=projection, lon_0=lon_0, resolution='c') # Mollweide projection, medium res
#         a.set_adjustable('box') # fix the aspect ratio, so resizing figure doesn't change box
#         m.drawcoastlines(color='#AAAAAA')
#         m.drawparallels(np.arange(-90,91,30),color='#AAAAAA')
#         m.drawmeridians(np.arange(-135,136,45),color='#AAAAAA')
#     return m # for coordinate conversions


#------------------------------------------------------------------------------
# Classes
#------------------------------------------------------------------------------
class Normalize(mpl.colors.Normalize):
    """
    Class for passing into kwarg norm==, when declaring new pcolor or 
    contourf objects. Creates new normalizations of colormap.
    
    TODO: FIX THIS! IDEALLY, COLORS ON COLORBAR HAVE SAME GRADIATION BUT 
    XTICKLABELS ARE WARPED; INSTEAD, COLORBAR COLORS SEEM TO STRETCH AND TICKS
    ARE STILL LINEAR. OR... MAYBE THAT'S WHAT WE WANT.
    """
    def __init__(self, exp=0, midpoint=None, vmin=None, vmax=None, clip=None):
        # User will use -10 to 10 scale; converted to value used in equation
        if abs(exp) > 10: raise ValueError('Warping scale must be between -10 and 10.')
        self.midpoint = midpoint
        self.exp = exp
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Function
        def warp(x, exp, exp_max=4):
            # Returns indices stretched so neutral/low values are sampled more heavily
            # ...will artifically use exp to signify stretching away from neutral vals,
            # ...or compressing toward neutral vals
            if exp > 0:
                invert = True
            else:
                invert, exp = False, -exp
            exp = exp*(exp_max/10)
            # ...apply function; approaches x=1 as a-->Inf, x=x as a-->0
            if invert: x = 1-x
            value =  (x-1+(np.exp(x)-x)**exp)/(np.e-1)**exp
            if invert: value = 1-value # flip on y-axis
            # ...and return
            return value
        # Initial stuff
        if self.midpoint is None:
            midpoint = self.vmin
        else:
            midpoint = self.midpoint
        # Get middle point in 0-1 coords, and value
        midpoint_scaled = (midpoint - self.vmin)/(self.vmax - self.vmin)
        value_scaled = (value - self.vmin)/(self.vmax - self.vmin)
        try:
            iter(value_scaled)
        except TypeError:
            value_scaled = np.arange(value_scaled)
        value_cmap = np.ma.empty(value_scaled.size)
        for i,v in enumerate(value_scaled):
            # ...get values, accounting for midpoints
            if v < 0: v = 0
            if v > 1: v = 1
            if v >= midpoint_scaled:
                block_width = 1 - midpoint_scaled
                value_cmap[i] = (midpoint_scaled + 
                        block_width*warp((v - midpoint_scaled)/block_width, self.exp)
                        )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled - 
                        block_width*warp((midpoint_scaled - v)/block_width, self.exp)
                        )
        return value_cmap

class Midpoint(Normalize):
    """
    Function for normalizing colormaps around midpoint.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

