#------------------------------------------------------------------------------
# All
#------------------------------------------------------------------------------
__all__ = [
    'Figure',
    'subplots', 'cmapload', # basic setup
    'Normalize', 'Formatter', 'FracFormatter', # custom colormap normalization
    ]

#------------------------------------------------------------------------------
# Imports, all
#------------------------------------------------------------------------------
from . import utils
import os # need for saving figures
import fractions # from standard lib
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt # will attach a suitable backend, make available the figure/axes modules
import mpl_toolkits.basemap as basemap
import numpy as np
import pandas as pd
# import cartopy.crs as ccrs # crs stands for "coordinate reference system", leading c is "cartopy"
# import cartopy.mpl as cmpl

#------------------------------------------------------------------------------
# Important class
#------------------------------------------------------------------------------
class Figure(mpl.figure.Figure):
    """
    Subclass of the mpl.figure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, *args, **kwargs):
        # mpl.figure.Figure.__init__(self, *args, **kwargs)
        # super(Figure, self).__init__(*args, **kwargs)
        super().__init__(*args, **kwargs) # python 3 only
        self.basemaps = [] # new atribute
        self.panelx = None # replace with axis object
        self.panely = None # replace with axis object

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
        TODO: 
            * Add options for labels spanning multiple subplots. Possibly forget some keyword groups;
            maybe just take 'color' and 'linewidth' for axis limits and fonts, with grid.
            * Add options for special legend handling; allow overwriting the label linewidth, etc.
            if we want to differentiate colors, for example.
            * Add options for datetime handling; note possible date axes objects are TimeStamp (pandas),
            np.datetime64, DateTimeIndex; can fix with fig.autofmt_xdate() or manually set options; uses
            ax.is_last_row() or ax.is_first_column(), so should look that up...problem is there is no 
            autofmt_ydate(), so really should implement my own version of this...
            * Add options/considerations for twinx/twiny; since this iterates through all axes children in the 
            figure object, should be valid. Should be **added to my.subplots**, actually, and then
            my.subplots can handle their placement in the array.
            * Write module for custom lat/lon annotations... should be able to be done
            with .annotate specifying lon/lat coordinates, I believe.
            * Override basemap methods so that I DON'T HAVE TO TO PADDING IN LONGITUDE, AND DON'T
            HAVE TO USE MESHGRID EVERYTIME, AND DON'T HAVE TO USE LATLON=TRUE...or can have it EMPLOY
            SEAMFIX METHOD, WHICH IS ALSO CAPABLE OF CIRCULARLY SHIFTING LONGITUDES

        NOTE: sharex and sharey only affect *ticklocations*/*ticklabels*/*tickformatters*/*limits*
        by making axes correspond to each other. Does *not* affect axis labels, unless you used
        my.subplots (which set them invisible); this function also makes them invisible.
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
        legend_kw = utils.dict(dict(framealpha=0.6, fancybox=False), **legend_kw)
        colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=False), **colorbar_kw)
        # ...misc
        abc_convert = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e', 5:'f', 6:'g', 7:'h', 8:'j', 9:'i', 10:'k', 11:'l', 12:'m'}

        # Colorbar/legend checks
        if (colorbar or colorbar_ax is not None) and mappable is None: # steal space from existing axis
            raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
        if legend_ax is not None and objects is None:
            raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')

        # Set properties
        for i,a in enumerate(self.axes): # .flat gives row-major indexing, just like matlab for plots
            # Special axes testing
            try: legend_ax_on = ((legend_ax % len(self.axes))==i)
            except TypeError: legend_ax_on = False
            try: colorbar_ax_on = ((colorbar_ax % len(self.axes))==i)
            except TypeError: colorbar_ax_on = False
            basemap_ax_find = np.where([a==basemap.ax for basemap in self.basemaps])[0]

            # Legend axes setup
            # ...some checks
            if legend_ax_on: # e.g. can input -1, or axes.Axes object
                legend_ax_kw = dict(legend_kw, handles=objects, loc='center')
                a.set_frame_on(False)
                a.xaxis.set_visible(False)
                a.yaxis.set_visible(False)
            elif legend:
                legend_ax_kw = legend_kw
            if legend or legend_ax_on:
                # ...draw legend
                leg = a.legend(**legend_ax_kw)
                # ...check for certain properties only
                mpl.pyplot.setp(leg.get_frame(), **patch_kw)
                # ...and apply text properties, identical to ticklabel properties
                mpl.pyplot.setp(leg.get_texts(), **ticklabels_kw) # color should be in there as font property
                if legend_ax_on:
                    continue

            # Colorbar axes setup
            # ...some checks
            if colorbar_ax_on:
                colorbar_ax_kw = dict(colorbar_kw, cax=a, use_gridspec=True) # uses space affored by this entire axis
            elif colorbar:
                colorbar_ax_kw = dict(colorbar_kw, ax=a) # steals space form axis 'a', uses make_axes() to build new one
            if colorbar or colorbar_ax_on:
                # ...draw colorbar
                cb = self.colorbar(mappable, **colorbar_ax_kw)
                # ...first, make edges/dividers consistent with axis edges
                if cb.dividers is not None: mpl.pyplot.setp(cb.dividers, **axes_kw) # change dividers
                mpl.pyplot.setp(cb.outline, **patch_kw)
                # ...also take care of pesky white lines between color segments
                cb.solids.set_edgecolor('face')
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
                if colorbar_ax_on:
                    continue

            # Title, and abc'ing, and suptitle
            if i==0 and suptitle is not None:
                self.suptitle(suptitle, **title_kw)
            if titles is not None:
                a.set_title(titles[i], **title_kw)
            elif title is not None:
                a.set_title(title, **title_kw)
            if abc:
                a.annotate(abc_convert[i] + ')', xy=(0,0), xytext=(0,1.03),
                        xycoords='axes fraction', textcoords='axes fraction',
                        verticalalignment='baseline', horizontalalignment='left',
                        **title_kw)

            # Projection axes setup
            # print(basemap_ax_find)
            if basemap_ax_find.size>0:
                # Coatlines
                m = self.basemaps[basemap_ax_find[0]]
                m.drawcoastlines(**coastlines_kw)
                # Parallels
                pl, md = m.drawparallels(utils.arange(-90,90,30)), m.drawmeridians(utils.arange(-180,180,60))
                for pi in pl.values():
                    mpl.pyplot.setp(pi, **gridlines_kw)
                for mi in md.values():
                    mpl.pyplot.setp(mi, **gridlines_kw)
                # Frame boundary (forgot to set that)
                # ...add stuff here
                continue # skip everything else

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
            if ylabels is not None:
                ylabel = ylabels[i]
            if xlabels is not None:
                xlabel = xlabels[i]
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
            
#------------------------------------------------------------------------------
# Definitions
#------------------------------------------------------------------------------
def subplots(array=None, nrows=1, ncols=1, sharex=False, sharey=False, # for sharing x, y axes
        aspect=None, height=None, width=5, # for controlling aspect ratio; default is control for width
        panelx=False, panely=False, panelx_width=0.2, panely_width=0.2, # add panels for legends, colorbars, etc.
        hspace=0.5, wspace=0.5, height_ratios=None, width_ratios=None, # options for gridspec (instead of gridpsec_kw)
            # spacing betwen axes, in inches
        left=0.5, bottom=0.5, right=0.2, top=0.3,
            # edges of plot for GridSpec object, in inches
        toolkit='basemap', projection=None, **basemap_kw): # for projections; can be 'basemap' or 'cartopy'
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want a smallish panel
    to plot data, just use width/height ratios. don't bother with kwarg dictionaries because
    it is really, really customized; we do use them in formatting function though.
    
    TODO:
        * if want **differing spacing** between certain subplots, must use 
        GridSpecFromSubplotSpec with its **own** hspace/wspace properties...requires
        some consideration, but don't worry until you come across it
        * figure out how to automate fixed aspect ratios for basemap subplots, so
        don't have to guess/redraw a bunch of times
        * figure size should be constrained by WIDTH/HEIGHT/ASPECT RATIO OF AXES, 
        with everything else provided; must pick two of these; default, is to pick
        width and aspect, but if height provided, width is ignored instead. no way
        right now to determine aspect by constraining width/height
        * Use GridSpecFromSubplotSpec to ENABLE DIFFERENG WSPACES/HSPACES IN A SUBPLOT GRID, FOR EXAMPLE, 
        FOR A COLORBAR LYING ON THE BOTTOM OF A PLOT; top/bottom/left/right should be zero here, but don't have
        to worry; function doesn't even take them as arguments...should be easy
    """
    # First things first, interpret rarray input and ratios
    # ...array stuff
    if array is None:
        array = np.reshape(np.arange(1,nrows*ncols+1), (nrows, ncols))
    array = np.array(array)
    nrows, ncols = array.shape[0], array.shape[1]
        # ...will mean NUMBER, NOT INCLUDING PANELS
    # ...stuff for getting shared axes
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # ...0 or -1 stands for 'empty'; could also use NaN
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    num_axes = len(axes_ids)
    # ...width/height ratios
    if width_ratios is None:
        width_ratios = np.ones(ncols)/ncols
    else:
        width_ratios = np.array(width_ratios)/sum(width_ratios)
    if height_ratios is None:
        height_ratios = np.ones(nrows)/nrows
    else:
        height_ratios = np.array(height_ratios)/sum(height_ratios)
    # ...for projections (default determine aspect ratio)
    if projection is not None and aspect is None:
        aspects = (yrange[:,1]-yrange[:,0])/(xrange[:,1]-xrange[:,0])
        if not np.all(aspects[0]==aspects): # should ammend to accept width/height ratios
            raise ValueError('If each subplot is map projection, their aspect ratios should match.')
        aspect = .5*aspects[0] # .6 is good bet
    
    # Fixed aspect, infer required figure size [try wspace=0.2 for gs with cols=4, axes gs[0,:2] and gs[0,2:] 
    # vs. gs with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ])
    # FORMULA: aspect = ((width - left - right - (ncol-1)*wspace)/ncol) / ((height - top - bottom - (nrow-1)*hspace)/nrow)
    # ...set aspect
    if aspect is None: # is still empty, that is
        aspect = 1
    # ...set figure aspect ratio, fixing height
    if height is not None:
        # height = width*(aspect*(ncols + (ncols-1)*hspace)) / (nrows + (nrows-1)*wspace)
        # axheight_ave_nopanel = (height - top - bottom - (nrows-1)*hspace - int(panelx)*(panelx_width + hspace))/nrows
        axheight_ave_nopanel = (height - top - bottom - (nrows-1)*hspace - int(panelx)*(panelx_width + bottom))/nrows # SUBPLOTSPEC
        if axheight_ave_nopanel<0: raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axwidth_ave_nopanel = (axheight_ave_nopanel*aspect)
        # width = axwidth_ave_nopanel*ncols + left + right + (ncols-1)*wspace + int(panely)*(panely_width + wspace)
        width = axwidth_ave_nopanel*ncols + left + right + (ncols-1)*wspace + int(panely)*(panely_width + right) # SUBPLOTSPEC
    # ...set figure aspect ratio, fixing width
    else:
        # width = height*(nrows + (nrows-1)*wspace) / (aspect*(ncols + (ncols-1)*hspace))
        # axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - int(panely)*(panely_width + wspace))/ncols # average in inches
        axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - int(panely)*(panely_width + right))/ncols # SUBPLOTSPEC
        if axwidth_ave_nopanel<0: raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axheight_ave_nopanel = (axwidth_ave_nopanel/aspect)
        # height = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + int(panelx)*(panelx_width + hspace)
        height = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + int(panelx)*(panelx_width + bottom) # SUBPLOTSPEC
    figsize = (width, height)
    # ...also, don't worry about left/bottom/top/right, because we will save figure
    # with tight bounding box anyway... or not

    # Convert hspace/wspace/edges to INCHES units; in the above formula, they had to all be inches, note
    # ...first, width
    axwidth_ave_panel = (axwidth_ave_nopanel*ncols + int(panely)*panely_width)/(ncols + int(panely))
    # wspace = wspace/axwidth_ave_panel
    wspace = wspace/axwidth_ave_nopanel # SUBPLOTSPEC
    # ...next, height
    axheight_ave_panel = (axheight_ave_nopanel*nrows + int(panelx)*panelx_width)/(nrows + int(panelx))
    # hspace = hspace/axheight_ave_panel
    hspace = hspace/axheight_ave_nopanel # SUBPLOTSPEC
    # ...next, width/height ratios and panel stuff
    if panely: # right panel
        # width_ratios = np.append(width_ratios, panelx_width/axwidth_ave_nopanel)
        # width_ratios = np.append(width_ratios, panelx_width/(ncols*axwidth_ave_nopanel))
        # width_ratios = width_ratios/width_ratios.sum()
        # array = np.concatenate((array, np.ones((nrows, 1))*(array.max()+1)), axis=1)
        wspace_outer = right/((ncols*axwidth_ave_nopanel + (ncols-1)*wspace + panely_width)/2) # SUBPLOTSPEC; the left "axis" is the main gridspec object including wspaces
        width_ratios_outer = np.array([ncols*axwidth_ave_nopanel + (ncols-1)*wspace, panely_width])
        width_ratios_outer = width_ratios_outer/width_ratios_outer.sum()
    else:
        wspace_outer, width_ratios_outer = None, None
    if panelx: # bottom panel
        # height_ratios = np.append(height_ratios, panely_width/axheight_ave_nopanel)
        # height_ratios = np.append(height_ratios, panely_width/(nrows*axheight_ave_nopanel))
        # height_ratios = height_ratios/height_ratios.sum()
        # column = np.ones((1, ncols+int(panely)))*(array.max()+1)
        # if panely: column[-1] = 0
        # array = np.concatenate((array, column), axis=0)
        hspace_outer = bottom/((nrows*axheight_ave_nopanel + (nrows-1)*hspace + panelx_width)/2) # SUBPLOTSPEC; the top "axis" is the main gridspec object including hspaces
        height_ratios_outer = np.array([nrows*axheight_ave_panel + (nrows-1)*hspace, panelx_width])
        height_ratios_outer = height_ratios_outer/height_ratios_outer.sum()
    else:
        hspace_outer, height_ratios_outer = None, None
    # ...and finally, edges
    left, right = left/width, 1-right/width
    bottom, top = bottom/height, 1-top/height

    # Set up figure and gridspec object, with custom figure class
    # fig = mpl.pyplot.figure(FigureClass=Figure, figsize=figsize)
    # gs = mpl.gridspec.GridSpec(nrows=nrows+int(panelx), ncols=ncols+int(panely), # add row/column if panel exists
    #         left=left, bottom=bottom, right=right, top=top,
    #         wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)
    
    # Or, set up with subplotspec
    fig = mpl.pyplot.figure(FigureClass=Figure, figsize=figsize)
    if panelx or panely:
        GS = mpl.gridspec.GridSpec(nrows=1+int(panelx), ncols=1+int(panely), left=left, bottom=bottom, right=right, top=top,
                wspace=wspace_outer, hspace=hspace_outer, width_ratios=width_ratios_outer, height_ratios=height_ratios_outer) # set wspace/hspace to match the top/bottom spaces
        gs = mpl.gridspec.GridSpecFromSubplotSpec(nrows=nrows, ncols=ncols, subplot_spec=GS[0,0],
                wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)
        print(wspace_outer, hspace_outer, width_ratios_outer, height_ratios_outer)
    else:
        gs = mpl.gridspec.GridSpec(nrows=nrows+int(panelx), ncols=ncols+int(panely), # add row/column if panel exists
                left=left, bottom=bottom, right=right, top=top,
                wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)
    
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
    # ...create panel axes (at end)
    if panely:
        fig.panely = fig.add_subplot(GS[0,1])
        # a = fig.add_subplot(gs[:,-1])
        # ax = np.append(ax, a)
    if panelx:
        fig.panelx = fig.add_subplot(GS[1,0])
        # a = fig.add_subplot(gs[-1,:])
        # ax = np.append(ax, a)
    # ...set basemap
    if projection is not None:
        m = np.array(num_axes*[None])
        for i in range(num_axes):
            m[i] = basemap.Basemap(projection=projection, ax=ax[i], **basemap_kw)
        fig.basemaps = m.tolist()
    # ...take out axes
    if ax.size==1:
        ax = ax[0]
    if projection is not None:
        if m.size==1:
            m = m[0]

    # Finally, return result
    print('Axes grid constructed.')
    if projection is None:
        return fig, ax
    else:
        return fig, ax, m

def cmapload(name, directory='~/cmaps'):
    """
    Load colormaps from RGB file using pandas.
    Call this function whenever startying up matplotlib.
    Note user can use ~; pandas should replace it with the $HOME directory.
    """
    # Directory
    has_neutral, has_cutoff = False, False
    if directory[-1]!='/':
        directory = directory + '/'
    
    # Load
    if type(name) is str:
        if name=='bw':
            cmap = np.tile(utils.arange(0, 1, 0.1)[np.newaxis,:], (1,3))
        # Colorbrewer
        elif name=='hclPurple':
            cmap = pd.read_table(directory + 'hclPurple.rgb', sep=',', header=None, skiprows=1).values
        elif name=='hclHot':
            cmap = pd.read_table(directory + 'hclHot.rgb', sep=',', header=None, skiprows=1).values
        elif name=='hclRed':
            cmap = pd.read_table(directory + 'hclRed.rgb', sep=',', header=None, skiprows=1).values
        elif name=='hclWet':
            cmap = pd.read_table(directory + 'hclWet.rgb', sep=',', header=None, skiprows=1).values
        elif name=='hclPretty':
            cmap = pd.read_table(directory + 'hclPretty.rgb', sep=',', header=None, skiprows=1).values
        elif name=='hclCool':
            cmap = pd.read_table(directory + 'hclCool.rgb', sep=',', header=None, skiprows=1).values
        # Rainbows
        elif name=='rainbow':
            cmap = pd.read_table(directory + 'WhiteBlueGreenYellowRed.rgb', sep='\s+', header=None).values
            cmap = cmap/255
        # Blue-reds: hard
        elif name=='bwo': # hot reddish-orange to blue, white in middle
            cmap = pd.read_table(directory + 'BlueWhiteOrangeRed.rgb', sep='\s+', header=None).values
            cmap = cmap/255 # has white in middle, blue to light blue to white, yellow, red
            has_neutral = True
        elif name=='br':
            cmap = pd.read_table(directory + 'BlRe.rgb', sep='\s+', header=None).values
            cmap = cmap/255 # red-blue very saturated, with sharp transition
            has_cutoff = True
        elif name=='brsoft':
            cmap = pd.read_table(directory + 'BlueRed.rgb', skiprows=1, sep='\s+', header=None).values
            cmap = cmap/255 # NEW
            has_cutoff = True
        elif name=='bo':
            cmap = pd.read_table(directory + 'BlueYellowRed.rgb', sep='\s+', header=None).values
            cmap = cmap/255 # blue to orange/red, with SHARP transition in middle
            has_cutoff = True
        elif name=='wbrw':
            cmap = pd.read_table(directory + 'WhBlReWh.rgb', skiprows=1, sep='\s+', header=None).values
            cmap = cmap/255 # NEW
            has_cutoff = True
        # Land
        elif name=='terrain':
           cmap = pd.read_table(directory + 'OceanLakeLandSnow.rgb', sep='\s+', header=None).values
           cmap = cmap[2:,:]/255 # ignore the ocean/lake colors
        elif name=='night':
            cmap = pd.read_table(directory + 'GMT_nighttime.rgb', skiprows=2, sep='\s+', header=None).values
            has_cutoff = True
        # Dry-wets
        elif name=='drywet':
            cmap = pd.read_table(directory + 'GMT_drywet.rgb', skiprows=2, sep='\s+', header=None).values
        elif name=='wet': # very few levels
            cmap = pd.read_table(directory + 'CBR_wet.rgb', skiprows=2, sep='\s+', header=None).values
            cmap = cmap/255
        elif name=='Wet': # very few levels
            cmap = pd.read_table(directory + 'precip_11lev.rgb', skiprows=2, sep='\s+', header=None).values
            cmap = cmap/255
        elif name=='DryWet': # brown white green blue
            cmap = pd.read_table(directory + 'precip_diff_12lev.rgb', skiprows=2, sep='\s+', header=None).values
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
    
    # Number of colors by default
    N = cmap.shape[0]

    # And create object; need ListedColormap here, not LinearSegmentedColormap (see documentation)
    # cmap_dict = utils.dict(red=cmap[:,0], green=cmap[:,1], blue=cmap[:,2])
    # if cmap.shape[1]==4: cmap_dict.alpha = cmap[:,3]
    # cmap_object = mpl.colors.LinearSegmentedColormap(name, cmap_dict, cmap.shape[0])
    # cmap_object = mpl.colors.ListedColormap(cmap, name=name, N=cmap.shape[0])
    
    # Or make colormap that linearly interpolates between each registered color
    cmap_object = mpl.colors.LinearSegmentedColormap.from_list(colors=cmap, name=name, N=256)

    # Register, and return object (user can call by name, or by supplying cmap=cmap_object
    # to a contouring/pcolor function)
    mpl.pyplot.register_cmap(cmap=cmap_object)
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
# Classes and Formatters
#------------------------------------------------------------------------------
class Basemap(basemap.Basemap):
    """
    TODO: WRITE ENTIRE SUBCLASS.
    """
    def __init__(self, *args, **kwargs):
        # basemap.Basemap.__init__(self, *args, **kwargs)
        # super(Basemap, self).__init__(*args, **kwargs)
        super().__init__(*args, **kwargs)
        if 'lon_0' in kwargs:
            self.lon_0 = lon_0 # by default, only lonmin/lonmax are written
    
    # Seam-fixing functions
    def wrapper(lon, lat, data, lon_0=-180):
        """
        Adjusts plot for basemap plotting, with seams.
        TODO: WRITE ENTIRE FUNCTION; should just use AVERAGE OF SUB-POLAR CELLS ==
        THE 90/-90 VALUE, and also WRAP AROUND
        """
        # Note data should be lat by lon for plotting commands
        # ...fix longitudes
        lon, data = geo.lonfix(lon, data.T, lon_0=lon_0)
        data = data.T
        # ...get seam values
        weights = (lon_0+360-lon[-1], lon_0-lon[0])
        data_seam = weights[0]*data[:,-1:] + weights[1]*data[:,:1]
        data = np.concatenate((data_seam, data, data_seam), axis=1)
        # ...fix poles (simple)
        data_sp = np.tile(data[-1:,:].mean(axis=1, keepdims=True), (1,data.shape[1]))
        data_np = np.tile(data[:1,:].mean(axis=1, keepdims=True), (1,data.shape[1]))
        data = np.concatenate((data_sp, data, data_np), axis=0)
        # ...and lons/lats
        lon = np.append(np.append(lon_0, lon), lon_0+360)
        lat = np.append(np.append(-90, lat), 90)
        Lon, Lat = np.meshgrid(lon, lat)
        return Lon, Lat, data

    # Re-define plotting comands
    def contourf(self, lon, lat, data, **kwargs):
        lon, lat, data = wrapper(lon, lat, data)
        # super(Basemap, self).contourf(lon, lat, data, lonlat=True, **kwargs)
        super().contourf(lon, lat, data, lonlat=True, **kwargs)

class Normalize(mpl.colors.Normalize):
    """
    Class for passing into kwarg norm==, when declaring new pcolor or 
    contourf objects. Creates new normalizations of colormap.
    
    TODO: FIX THIS! IDEALLY, COLORS ON COLORBAR HAVE SAME GRADIATION BUT 
    XTICKLABELS ARE WARPED; INSTEAD, COLORBAR COLORS SEEM TO STRETCH AND TICKS
    ARE STILL LINEAR. OR... MAYBE THAT'S WHAT WE WANT.
    """
    def __init__(self, exp=0, extend='neither', midpoint=None, vmin=None, vmax=None, clip=None):
        # User will use -10 to 10 scale; converted to value used in equation
        if abs(exp) > 10: raise ValueError('Warping scale must be between -10 and 10.')
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        self.exp = exp
        self.extend = extend
        # mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

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
        if self.extend=='both' or self.extend=='max':
            value_cmap[value_cmap>1] = 1
        if self.extend=='both' or self.extend=='min':
            value_cmap[value_cmap<0] = 0
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

def Formatter(sigfig=3):
    """
    Format just like a number.
    """
    # Format definition
    def f(n, loc, sigfig=sigfig):
        formatstr = '%%.%.ff' % (sigfig,)
        string = formatstr % (n,)
        return string.rstrip('0').rstrip('.')
    # And create object
    return mpl.ticker.FuncFormatter(f)

def FracFormatter(fact=np.pi, symbol=r'\pi'):
    """
    Formatter for fractions along axis.
    """
    # Start with fraction definition
    def frac(n, loc, fact=fact, symbol=symbol): # must accept location argument
        f = fractions.Fraction(n/fact).limit_denominator()
        if n==0: # zero
            return '0'
        elif f.denominator==1: # denominator is one
            if f.numerator==1:
                return r'$%s$' % (symbol,)
            elif f.numerator==-1:
                return r'${-}%s$' % (symbol,)
            else:
                return r'$%d%s$' % (f.numerator, symbol)
        elif f.numerator==1: # numerator is +/-1
            return r'$%s/%d$' % (symbol, f.denominator)
        elif f.numerator==-1:
            return r'${-}%s/%d$' % (symbol, f.denominator)
        else: # otherwise
            return r'$%d%s/%d$' % (f.numerator, symbol, f.denominator)
    # And create FuncFormatter class
    return mpl.ticker.FuncFormatter(frac)

