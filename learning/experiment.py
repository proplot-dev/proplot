# Experiment with subclassing the figure class
# ...should probably just make a simple class that is completely identical to myformat, 
# except it iterates over the axes in axes.
import matplotlib as mpl
import numpy as np
import utils

#------------------------------------------------------------------------------
# Subclasses for custom, simple formatting options
#------------------------------------------------------------------------------
class Figure(mpl.figure.Figure):
    """
    Subclass of the mpl.figure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, *args, **kwargs):
        mpl.figure.Figure.__init__(self, *args, **kwargs)

    def format(self, filename=None, desktop=True, # file name, and whether to put in desktop folder
            m=None, # basemap options; the basemap object
            legend=False, objects=None, labels=None, legend_ax=None, legend_kw=dict(), # legend options; if declared and have no objects with 'label', does nothing
            colorbar=False, mappable=None, colorbar_ax=None, colorbar_kw=dict(), # colorbar options
            xgrid=True, ygrid=True, xgridminor=True, ygridminor=True, # grid options (add to this)
            xtickminor=False, ytickminor=False, # turning on/off minor ticks
            xtickpos='bottom', ytickpos='left', xopposite=True, yopposite=True, # location of ticks, which spines to draw
            xlim=None, ylim=None, xscale='linear', yscale='linear', xreverse=False, yreverse=False, # special properties
            xdates=False, ydates=False, # whether to format axis labels as long datetime strings
            xlabel=None, ylabel=None, clabel=None, title=None, suptitle=None, abc=False, # axis/colorbar labels, titles
            xlabels=None, ylabels=None, titles=None, # options for multiple labeling
            xlocator=None, xminorlocator=None, clocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
            xformatter=None, cformatter=None, yformatter=None, # formatter
            axes_kw = dict(), # for color/linewidth
            ctick_kw=dict(), xtick_kw=dict(), ytick_kw=dict(), # all should be special tick properties, using set_tick_params; also takes axis_kw color
            ticklabels_kw=dict(), label_kw=dict(), title_kw=dict(), suptitle_kw=dict(), abc_kw=dict(), # all are passed to FontProperties class; also takes axis_kw color
            grid_kw=dict(), gridminor_kw=dict(), gridlines_kw=dict(), coastlines_kw=dict(), # for axis properties, all are line/patch-type properties
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

        NOTE: title, xlabel, and ylabel can all be multi-dimensional
        """
        # Override defaults for legend/colorbar, labels, ticklabels, title, grid/gridminor
        # ...axes spines, ticks (set with set_tick_params)
        axes_kw = utils.dict(dict(color='k', linewidth=1), **axes_kw)
        xtick_kw = utils.dict(dict(length=5, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **xtick_kw)
        ytick_kw = utils.dict(dict(length=3, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **ytick_kw)
        ctick_kw = utils.dict(dict(length=5, direction='out', width=axes_kw.linewidth, color=axes_kw.color), **ctick_kw)
        patch_kw = utils.dict(linewidth=axes_kw.linewidth, edgecolor=axes_kw.color)
        # ...grid lines (set on declaration)
        grid_kw = utils.dict(dict(linestyle='-', linewidth=0.8, color=axes_kw.color, alpha=0.1), **grid_kw)
        gridminor_kw = utils.dict(dict(linestyle=':', linewidth=0.5, color=axes_kw.color, alpha=0.05), **gridminor_kw)
        xtickminor_kw = utils.dict(xtick_kw, length=xtick_kw.length/2) # dependent on above
        ytickminor_kw = utils.dict(ytick_kw, length=ytick_kw.length/2)
        # ...geographic (set on declaration)
        coastlines_kw = utils.dict(dict(linewidth=axes_kw.linewidth, color=axes_kw.color), **coastlines_kw)
        gridlines_kw = utils.dict(dict(linestyle=':', linewidth=axes_kw.linewidth, color=axes_kw.color, alpha=0.2), **gridlines_kw)
        # ...text (set with plt.setp, set_label_text, or on creation)
        ticklabels_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **ticklabels_kw)
        label_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **label_kw)
        title_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **title_kw)
        abc_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **abc_kw)
        # ...special (set on creation)
        legend_kw = utils.dict(dict(framealpha=1, fancybox=False), **legend_kw)
        colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=False), **colorbar_kw)
        # ...misc
        abc_convert = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e', 5:'f', 6:'g', 7:'h', 8:'j', 9:'i', 10:'k', 11:'l', 12:'m'}

        # Colorbar/legend checks
        if colorbar and mappable is None: # steal space from existing axis
            raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
        if legend and objects is None and legend_ax is not None:
            raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')

        # Set properties
        for i,a in enumerate(self.axes): # .flat gives row-major indexing, just like matlab for plots
            print('Axes %d setup...' % i)
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
                mpl.pyplot.setp(leg.get_frame(), **patch_kw)
                # ...and apply text properties, identical to ticklabel properties
                mpl.pyplot.setp(leg.get_texts(), **ticklabels_kw) # color should be in there as font property

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
                cb = self.colorbar(mappable, **special_kw, **colorbar_kw)
                # ...first, make edges/dividers consistent with axis edges
                if cb.dividers is not None: mpl.pyplot.setp(cb.dividers, **axes_kw) # change dividers
                mpl.pyplot.setp(cb.outline, **patch_kw)
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
                mpl.pyplot.setp(axis.get_ticklabels(), **ticklabels_kw)
                axis.set_tick_params(which='both', **ctick_kw)
                if clabel is not None: axis.set_label_text(clabel, **label_kw)
                # ...how can i do this
                # ...problem is axis has vmin, vmax
                # ...and next, tick stuff (hard to deal with lines separating colors, vs. tick marks)

            # Projection axes setup
            if m is not None and not skip:
                skip = True
                m[i].drawcoastlines(**coastlines_kw)
                pl, md = m[i].drawparallels(utils.arange(-90,90,30)), m[i].drawmeridians(utils.arange(-180,180,60))
                for pi in pl.values():
                    mpl.pyplot.setp(pi, **gridlines_kw)
                for mi in md.values():
                    mpl.pyplot.setp(mi, **gridlines_kw)
                # ...set things up

            # Skip, if axes was either meant to be filled with legend or colorbar, or was projection
            if skip:
                continue

            # Deal with multiple titles/labels
            if titles is not None:
                title = titles[i]
            if ylabels is not None:
                ylabel = ylabels[i]
            if xlabels is not None:
                xlabel = xlabels[i]

            # Axes properties, that must be set in unison or have special properties
            # ...spine line properties, tick line properties (should be identical between axes)
            mpl.pyplot.setp(tuple(a.spines.values()), **axes_kw) # axis properties; universal
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
                mpl.pyplot.setp(axis.get_ticklabels(), **ticklabels_kw)
                # ...and finally, label
                if label is not None: axis.set_label_text(label)
                mpl.pyplot.setp(axis.get_label(), **label_kw)
                # ...and tick properties
                axis.set_tick_params(which='major', **tick_kw)
                axis.set_tick_params(which='minor', **tickminor_kw) # have length

            # Title, and abc'ing, and suptitle
            if i==0 and suptitle is not None:
                self.suptitle(suptitle, **title_kw)
            if title is not None:
                a.set_title(title, **title_kw)
            if abc:
                a.annotate(abc_convert[i] + ')', xy=(0,0), xytext=(0,1.03),
                        xycoords='axes fraction', textcoords='axes fraction',
                        verticalalignment='baseline', horizontalalignment='left',
                        **title_kw)

        # And save, optionally (into figures folder)
        if filename is not None:
            if '.' not in filename:
                filename = filename + '.pdf'
            if desktop:
                filename = os.environ['HOME'] + '/Desktop/Figures/' + filename
            print('Saving figure, filename: %s' % filename)
            self.savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.05, dpi=300,
                    transparent=True, frameon=False)
                # pad_inches is padding around bbox, and frameon makes figure window transparent
                # even though saving as pdf, think dpi will affect any embedded raster objects

def subplots(array=None, nrows=1, ncols=1, sharex=False, sharey=False, # for sharing x, y axes
        fixedaspect=False, aspect=None, height=None, width=7, # for controlling aspect ratio; default is control for width
        panelx=False, panely=False, panelx_ratio=0.05, panely_ratio=0.05, # add panels for legends, colorbars, etc.
        hspace=0.2, wspace=0.2, height_ratios=None, width_ratios=None, # options for gridspec (instead of gridpsec_kw)
        figsize=None, # options for plt.figure (right now, just figure size)
        projection=None, **basemap_kw): # for basemap
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want a smallish panel
    to plot data, just use width/height ratios. don't bother with kwarg dictionaries because
    it is really, really customized; we do use them in formatting function though.
    """
    # First things first, interpret rarray input and ratios
    #--------------------------------------------------------------------------
    # Array stuff
    #--------------------------------------------------------------------------
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
    # ...for projections
    if projection is not None:
        fixedaspect = True
        aspect = 0.8
    
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
    fig = mpl.pyplot.figure(FigureClass=Figure, figsize=figsize)
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
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])])
    if sharey:
        for i in ygroups_share:
            if ax[i] is None:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])])

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
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], sharex=sharex_ax, sharey=sharey_ax)
            # ...hide tick labels (not default behavior for manual sharex, sharey use)
            if sharex_ax is not None:
                mpl.pyplot.setp(ax[i].xaxis.get_ticklabels(), visible=False) 
                ax[i].xaxis.get_label().set_visible(False)
                # setp acts for iterable set of matplotlib objects; this is not black magic
                # these hold, because setting ticklabels or labels with higher-level commands
                # just changes the **text**, not underlying properties... pretty neat!
            if sharey_ax is not None:
                mpl.pyplot.setp(ax[i].yaxis.get_ticklabels(), visible=False)
                ax[i].yaxis.get_label().set_visible(False)
                
    # Extra stuff
    # ...set basemap
    if projection is not None:
        m = np.array(num_axes*[None])
        for i in range(num_axes):
            m[i] = basemap.Basemap(projection=projection, ax=ax[i], **basemap_kw)
    # ...create panel axes (at end)
    if panely:
        a = fig.add_subplot(gs[:,-1])
        ax = np.append(ax, a)
    if panelx:
        a = fig.add_subplot(gs[-1,:])
        ax = np.append(ax, a)

    # Finally, return result
    print('Axes grid constructed.')
    if projection is None:
        return fig, ax
    else:
        return fig, ax, m

# class Axes(mpl.axes.Axes):
#     """
#     Sublcass of the mpl.axes.Axes class, with some extra controls on legend/
#     colorbar options.
#     """
#     def __init__(self, *args, **kwargs):
#         super(MyAxes, self).__init__(*args, **kwargs)
#
#     # Add my axes
#     @staticmethod
#     def from_ax(ax=None, fig=None, *args, **kwargs):
#         if ax is None:
#             if fig is None: fig = figure(facecolor='w', edgecolor='w')
#             ax = MyAxes(fig, 1, 1, 1, *args, **kwargs)
#             fig.add_axes(ax)
#                 # ...and have to figure out how to add_subplot here
#         return ax
#
#     # Add series of functions for basemap plotting?
#     pass
#
# class Figure(mpl.figure.Figure):
#     """
#     Subclass of the mpl.figure.Figure class, with lots of special formatting
#     options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
#     in my subplots function.
#     """
#     def __init__(self, *args, **kwargs):
#         mpl.figure.Figure.__init__(self, *args, **kwargs)
#
#     def legend_format(self, ax, sample_kw=dict(), **kwargs):
#         """
#         Formats and declares legends; allow overwriting of the sample objects.
#         """
#         # Legend axes setup
#         # ...fist, make object
#         leg = ax.legend(**kwargs) # ...these are the legend kwargs
#         # ...check for certain properties only
#         plt.setp(leg.get_frame(), **patch_kw)
#         # ...and apply text properties, identical to ticklabel properties
#         plt.setp(leg.get_texts(), **text_kw) # color should be in there as font property
#         return
#
#     def colorbar_format(**kwargs):
#         """
#         Formats and declares colorbars.
#         """
#         # ...create object
#         cb = self.colorbar(mappable, **kwargs)
#         # ...first, make edges/dividers consistent with axis edges
#         if cb.dividers is not None: plt.setp(cb.dividers, **axes_kw) # change dividers
#         plt.setp(cb.outline, **patch_kw)
#         # ...next, ticks formatters/locators (different from axis_setup version)
#         if isinstance(clocator, mpl.ticker.Locator):
#             cb.locator = clocator
#         elif clocator is not None:
#             try: iter(clocator)
#             except TypeError: clocator = utils.arange(cb.vmin, cb.vmax, clocator)
#             cb.locator = mpl.ticker.FixedLocator(clocator)
#         if isinstance(cformatter, mpl.ticker.Formatter):
#             cb.formatter = cformatter
#         elif cformatter is not None:
#             if isinstance(cformatter, str):
#                 cb.formatter = mpl.ticker.FormatStrFormatter(cformatter) # %-style, numbers
#             else:
#                 cb.formatter = mpl.ticker.FixedFormatter(cformatter)
#         # ...and update (use this instead of update_bruteforce)
#         cb.update_ticks() # updates formatters/locators; colorbar uses weird combo of its own patch objects, and underlying axis
#         # ...finally, the ticks/ticklabels/labels basic properties
#         if cb.orientation=='horizontal':
#             axis = cb.ax.xaxis
#         else:
#             axis = cb.ax.yaxis
#         plt.setp(axis.get_ticklabels(), **ticklabels_kw)
#         axis.set_tick_params(which='both', **ctick_kw)
#         if clabel is not None: axis.set_label_text(clabel, **label_kw)
#         # ...how can i do this
#         # ...problem is axis has vmin, vmax
#         # ...and next, tick stuff (hard to deal with lines separating colors, vs. tick marks)
#         return
#
#     def format(self, filename=None, desktop=True, # file name, and whether to put in desktop folder
#             m=None, # basemap options; the basemap object
#             legend=False, objects=None, labels=None, legend_ax=None, legend_kw=dict(), # legend options; if declared and have no objects with 'label', does nothing
#             colorbar=False, mappable=None, colorbar_ax=None, colorbar_kw=dict(), # colorbar options
#             xgrid=True, ygrid=True, xgridminor=True, ygridminor=True, # grid options (add to this)
#             xtickminor=False, ytickminor=False, # turning on/off minor ticks
#             xtickpos='bottom', ytickpos='left', xopposite=True, yopposite=True, # location of ticks, which spines to draw
#             xlim=None, ylim=None, xscale='linear', yscale='linear', xreverse=False, yreverse=False, # special properties
#             xdates=False, ydates=False, # whether to format axis labels as long datetime strings
#             xlabel=None, ylabel=None, clabel=None, title=None, suptitle=None, abc=False, # axis/colorbar labels, titles
#             xlabels=None, ylabels=None, titles=None, # options for multiple labeling
#             xlocator=None, xminorlocator=None, clocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
#             xformatter=None, cformatter=None, yformatter=None, # formatter
#             axes_kw = dict(), # for color/linewidth
#             ctick_kw=dict(), xtick_kw=dict(), ytick_kw=dict(), # all should be special tick properties, using set_tick_params; also takes axis_kw color
#             ticklabels_kw=dict(), label_kw=dict(), title_kw=dict(), suptitle_kw=dict(), abc_kw=dict(), # all are passed to FontProperties class; also takes axis_kw color
#             grid_kw=dict(), gridminor_kw=dict(), gridlines_kw=dict(), coastlines_kw=dict(), # for axis properties, all are line/patch-type properties
#             ):
#         """
#         TODO: Add options for labels spanning multiple subplots.
#         Possibly forget some keyword groups; maybe just take 'color' and 'linewidth' for
#         axis limits and fonts, with grid.
#
#         TODO: Write module for custom lat/lon annotations... should be able to be done
#         with .annotate specifying lon/lat coordinates, I believe.
#
#         NOTE: possible date axes objects are TimeStamp (pandas), np.datetime64, DateTimeIndex;
#         can fix with fig.autofmt_xdate() or manually set options; uses ax.is_last_row()
#         or ax.is_first_column(), so should look that up...problem is there is no
#         autofmt_ydate(), so really should implement my own version of this...
#
#         NOTE: sharex and sharey only affect *ticklocations*/*ticklabels*/*tickformatters*/*limits*
#         by making axes correspond to each other. Does *not* affect axis labels, unless you used
#         my.subplots (which set them invisible); this function also makes them invisible.
#
#         NOTE: title, xlabel, and ylabel can all be multi-dimensional
#         """
#         # Override defaults for legend/colorbar, labels, ticklabels, title, grid/gridminor
#         # ...axes spines, ticks (set with set_tick_params)
#         axes_kw = utils.dict(dict(color='k', linewidth=1), **axes_kw)
#         xtick_kw = utils.dict(dict(length=5, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **xtick_kw)
#         ytick_kw = utils.dict(dict(length=3, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **ytick_kw)
#         ctick_kw = utils.dict(dict(length=5, direction='out', width=axes_kw.linewidth, color=axes_kw.color), **ctick_kw)
#         patch_kw = utils.dict(linewidth=axes_kw.linewidth, edgecolor=axes_kw.color)
#         # ...grid lines (set on declaration)
#         grid_kw = utils.dict(dict(linestyle='-', linewidth=0.8, color=axes_kw.color, alpha=0.1), **grid_kw)
#         gridminor_kw = utils.dict(dict(linestyle=':', linewidth=0.5, color=axes_kw.color, alpha=0.05), **gridminor_kw)
#         xtickminor_kw = utils.dict(xtick_kw, length=xtick_kw.length/2) # dependent on above
#         ytickminor_kw = utils.dict(ytick_kw, length=ytick_kw.length/2)
#         # ...geographic (set on declaration)
#         coastlines_kw = utils.dict(dict(linewidth=axes_kw.linewidth, color=axes_kw.color), **coastlines_kw)
#         gridlines_kw = utils.dict(dict(linestyle=':', linewidth=axes_kw.linewidth, color=axes_kw.color, alpha=0.2), **gridlines_kw)
#         # ...text (set with plt.setp, set_label_text, or on creation)
#         ticklabels_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **ticklabels_kw)
#         label_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **label_kw)
#         title_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **title_kw)
#         abc_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **abc_kw)
#         # ...special (set on creation)
#         legend_kw = utils.dict(dict(framealpha=1, fancybox=False), **legend_kw)
#         colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=False), **colorbar_kw)
#         # ...misc
#         abc_convert = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e', 5:'f', 6:'g', 7:'h', 8:'j', 9:'i', 10:'k', 11:'l', 12:'m'}
#
#         # Colorbar/legend checks
#         if colorbar and mappable is None: # steal space from existing axis
#             raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
#         if legend and objects is None and legend_ax is not None:
#             raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')
#
#         for i,a in enumerate(self.legend_axes):
#         # Format dedicated colorbar/legend axes
#         if self.colorbar_axes is not None:
#             self.colorbar_format(cax=self.colorbar_axes, use_gridspec=True, **colorbar_kw)
#         if self.legend_axes is not None:
#             self.legend_axes.set_frame_on(False)
#             self.legend_axes.xaxis.set_visible(False)
#             self.legend_axes.yaxis.set_visible(False)
#             self.legend_format(ax, **legend_kw)
#
#         for i,m in enumerate(self.basemaps):
#             # Projection axes setup
#             skip = True
#             m[i].drawcoastlines(**coastlines_kw)
#             pl, md = m[i].drawparallels(utils.arange(-90,90,30)), m[i].drawmeridians(utils.arange(-180,180,60))
#             for pi in pl.values():
#                 plt.setp(pi, **gridlines_kw)
#             for mi in md.values():
#                 plt.setp(mi, **gridlines_kw)
#             # ...set things up
#
#         # Set properties
#         for i,a in enumerate(self.axes): # .flat gives row-major indexing, just like matlab for plots
#             # Legend/colorbar for particular axes
#             if a.legend:
#                 legend_format(a)
#
#             # Deal with multiple titles/labels
#             if titles is not None:
#                 title = titles[i]
#             if ylabels is not None:
#                 ylabel = ylabels[i]
#             if xlabels is not None:
#                 xlabel = xlabels[i]
#
#             # Axes properties, that must be set in unison or have special properties
#             # ...spine line properties, tick line properties (should be identical between axes)
#             plt.setp(tuple(a.spines.values()), **axes_kw) # axis properties; universal
#             # ...tick alliegance
#             for axis, tickpos in zip((a.xaxis, a.yaxis), (xtickpos, ytickpos)):
#                 if tickpos is not None:
#                     axis.set_ticks_position(tickpos) # both/left/right
#             # ...spines visibility
#             for axis, sides, opposite in zip((a.xaxis, a.yaxis), (('bottom','top'), ('left','right')), (xopposite,yopposite)):
#                 for side in sides:
#                     if axis.get_ticks_position() in ('both',side) or opposite:
#                         a.spines[side].set_visible(True) # view top if ticks there, or turned it on
#                     elif not opposite:
#                         a.spines[side].set_visible(False)
#             # ...axis scaling
#             a.set_xscale(xscale)
#             a.set_yscale(yscale)
#             # ...axis limits
#             if xlim is None: xlim = a.get_xlim()
#             if ylim is None: ylim = a.get_ylim()
#             if xreverse: xlim = xlim[::-1]
#             if yreverse: ylim = ylim[::-1]
#             a.set_xlim(xlim)
#             a.set_ylim(ylim)
#             # ...and grid stuff
#             if xgrid:
#                 a.xaxis.grid(which='major', **grid_kw)
#             if ygrid:
#                 a.yaxis.grid(which='major', **grid_kw)
#             if xgridminor and xtickminor:
#                 a.xaxis.grid(which='minor', **gridminor_kw)
#             if ygridminor and ytickminor:
#                 a.yaxis.grid(which='minor', **gridminor_kw)
#
#             # X/Y Axis properties
#             for axis, label, dates, gridminor, tickminor, tickminorlocator, tickminor_kw, grid, ticklocator, tickformatter, tick_kw in zip(
#                     (a.xaxis, a.yaxis), (xlabel, ylabel), (xdates, ydates), # other stuff
#                     (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), (xtickminor_kw, ytickminor_kw), # minor ticks
#                     (xgrid, ygrid), (xlocator, ylocator), (xformatter, yformatter), (xtick_kw, ytick_kw) # major ticks
#                     ):
#                 # ...first, major tick locators (should not affect limits)
#                 lim = axis.get_view_interval() # to be used, automatically
#                 if isinstance(ticklocator, mpl.ticker.Locator):
#                     axis.set_major_locator(ticklocator)
#                 elif ticklocator is not None:
#                     try: iter(ticklocator)
#                     except TypeError: ticklocator = utils.arange(lim[0], lim[1], ticklocator)
#                     # axis.set_ticks(ticklocator)
#                     axis.set_major_locator(mpl.ticker.FixedLocator(ticklocator))
#                 # ...next, minor tick locator
#                 if not tickminor:
#                     axis.set_minor_locator(mpl.ticker.NullLocator())
#                 elif isinstance(tickminorlocator, mpl.ticker.Locator):
#                     axis.set_minor_locator(tickminorlocator)
#                 elif tickminorlocator is not None:
#                     try: iter(tickminorlocator)
#                     except TypeError: tickminorlocator = utils.arange(lim[0], lim[1], tickminorlocator)
#                     axis.set_minor_locator(mpl.ticker.FixedLocator(tickminorlocator))
#                 # ...next, tick formatters (enforce Null, always, for minor ticks)
#                 axis.set_minor_formatter(mpl.ticker.NullFormatter())
#                 if isinstance(tickformatter, mpl.ticker.Formatter): # use isinstance, because r"stuff" or f"stuff" **do not work**
#                     axis.set_major_formatter(tickformatter)
#                 elif tickformatter is not None:
#                     if isinstance(tickformatter, str):
#                         if dates: axis.set_major_formatter(mpl.dates.DateFormatter(tickformatter)) # %-style, dates
#                         else: axis.set_major_formatter(mpl.ticker.FormatStrFormatter(tickformatter)) # %-style, numbers
#                     else:
#                         # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
#                         axis.set_major_formatter(mpl.ticker.FixedFormatter(tickformatter)) # list of strings
#                 plt.setp(axis.get_ticklabels(), **ticklabels_kw)
#                 # ...and finally, label
#                 if label is not None: axis.set_label_text(label)
#                 plt.setp(axis.get_label(), **label_kw)
#                 # ...and tick properties
#                 axis.set_tick_params(which='major', **tick_kw)
#                 axis.set_tick_params(which='minor', **tickminor_kw) # have length
#
#             # Title, and abc'ing, and suptitle
#             if i==0 and suptitle is not None:
#                 a.get_figure().suptitle(suptitle, **title_kw)
#             if title is not None:
#                 a.set_title(title, **title_kw)
#             if abc:
#                 a.annotate(abc_convert[i] + ')', xy=(0,0), xytext=(0,1.03),
#                         xycoords='axes fraction', textcoords='axes fraction',
#                         verticalalignment='baseline', horizontalalignment='left',
#                         **title_kw)
#
#     # And save, optionally (into figures folder)
#     if filename is not None:
#         if '.' not in filename:
#             filename = filename + '.pdf'
#         if desktop:
#             filename = os.environ['HOME'] + '/Desktop/Figures/' + filename
#         print('Saving figure, filename: %s' % filename)
#         self.savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.05, dpi=300,
#                 transparent=True, frameon=False)
#             # pad_inches is padding around bbox, and frameon makes figure window transparent
#             # even though saving as pdf, think dpi will affect any embedded raster objects
#
