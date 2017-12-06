class Figure(mfigure.Figure):
    """
    Subclass of the mfigure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, *args, **kwargs):
        # Initialize
        # mfigure.Figure.__init__(self, *args, **kwargs)
        # super(Figure, self).__init__(*args, **kwargs)
        super().__init__(*args, **kwargs) # python 3 only
        # Advisories
        # self.axes_formatted = [] # nevermind
        # Special axes
        self.axes_main = []
        self.basemaps = [] # new atribute
        self.bottomcolorbar = None # overall colorbar; replace with axis object
        self.rightcolorbar = None
        self.bottomlegend = None # overall legend; replace with axis object
        self.title = None
    
    def save(self, filename, desktop=False):
        """
        Saving figures.
        """
        if '.' not in filename:
            filename = filename + '.pdf'
        if desktop: 
            filename = os.environ['HOME'] + '/Desktop/Figures/' + filename
        print('Saving figure, filename: %s' % filename)
        self.savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.05, dpi=300,
                transparent=True, frameon=False)
            # pad_inches is padding around bbox, and frameon makes figure window transparent
            # even though saving as pdf, think dpi will affect any embedded raster objects

    def format(self, # optionally specify the only axes to format?
            filename=None, desktop=True, # file name, and whether to put in desktop folder
            m=None, # basemap options; the basemap object
            legends=None, objects=None, labels=None, legend_ax=None, legend_kw=dict(), # legend options; if declared and have no objects with 'label', does nothing
            mappable=None, colorbar_kw=dict(), dividers_kw=dict(), # colorbar options; also colorbar dividers
            # colorbar=False, colorbar_ax=None,
            cgrid=False, xgrid=True, ygrid=True, xgridminor=True, ygridminor=True, # grid options (add to this); cgrid is the color dividers
            xtickminor=False, ytickminor=False, # turning on/off minor ticks
            xspineloc='both', yspineloc='both', xtickloc='bottom', ytickloc='top', # deals with spine options
            xlim=None, ylim=None, xscale='linear', yscale='linear', xreverse=False, yreverse=False, # special properties
            xdates=False, ydates=False, # whether to format axis labels as long datetime strings
            xlabel=None, ylabel=None, clabel=None, title=None, suptitle=None, abc=False, # axis/colorbar labels, titles
            xlabels=None, ylabels=None, titles=None, # options for multiple labeling
            xlocator=None, xminorlocator=None, clocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
            xformatter=None, cformatter=None, yformatter=None, # formatter
            # xlocators=None, xminorlocators=None, clocators=None, ylocators=None, yminorlocators=None, # optional differences
            # xformatters=None, cformatters=None, yformatters=None, # optional differences
            color='k', linewidth=0.7, transparent=True, # axes options
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
        xtick_kw = utils.dict(dict(length=5, direction='in', width=linewidth, color=color), **xtick_kw)
        ytick_kw = utils.dict(dict(length=3, direction='in', width=linewidth, color=color), **ytick_kw)
        ctick_kw = utils.dict(dict(length=5, direction='out', width=linewidth, color=color), **ctick_kw)
        spines_kw = utils.dict(linewidth=linewidth, color=color) # spine keyword
        patch_kw = utils.dict(linewidth=linewidth, edgecolor=color)
            # ...for outlines of colorbar, legend, cartopy frame, and cartopy features
        # ...grid lines (set on declaration)
        grid_kw = utils.dict(dict(linestyle='-', linewidth=0.8, color=color, alpha=0.1), **grid_kw)
        gridminor_kw = utils.dict(dict(linestyle=':', linewidth=0.5, color=color, alpha=0.05), **gridminor_kw)
        xtickminor_kw = utils.dict(xtick_kw, length=xtick_kw.length/2) # dependent on above
        ytickminor_kw = utils.dict(ytick_kw, length=ytick_kw.length/2)
        # ...geographic (set on declaration)
        coastlines_kw = utils.dict(dict(linewidth=linewidth, color=color), **coastlines_kw)
        gridlines_kw = utils.dict(dict(linestyle=':', linewidth=linewidth, color=color, alpha=0.2), **gridlines_kw)
        # ...text (set with plt.setp, set_label_text, or on creation)
        ticklabels_kw = utils.dict(dict(size=8, weight='normal', color=color), **ticklabels_kw)
        label_kw = utils.dict(dict(size=8, weight='normal', color=color), **label_kw)
        title_kw = utils.dict(dict(size=10, weight='normal', color=color), **title_kw)
        # title_kw = utils.dict(dict(size=10, weight='normal', color=color, loc='center'), **title_kw)
        abc_kw = utils.dict(dict(size=10, weight='bold', color=color), **abc_kw)
        # ...special (set on creation)
        legend_kw = utils.dict(dict(framealpha=0.6, fancybox=False), **legend_kw)
        colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=cgrid), **colorbar_kw) # use cgrid, for ease
        dividers_kw = utils.dict(dict(color=color, linewidth=linewidth), **dividers_kw)
        # ...misc
        abc_convert = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e', 5:'f', 6:'g', 7:'h', 8:'j', 9:'i', 10:'k', 11:'l', 12:'m'}
        # # ...axes spines, ticks (set with set_tick_params)
        # axes_kw = utils.dict(dict(color='k', linewidth=1), **axes_kw)
        # xtick_kw = utils.dict(dict(length=5, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **xtick_kw)
        # ytick_kw = utils.dict(dict(length=3, direction='in', width=axes_kw.linewidth, color=axes_kw.color), **ytick_kw)
        # ctick_kw = utils.dict(dict(length=5, direction='out', width=axes_kw.linewidth, color=axes_kw.color), **ctick_kw)
        # patch_kw = utils.dict(linewidth=axes_kw.linewidth, edgecolor=axes_kw.color)
        # # ...grid lines (set on declaration)
        # grid_kw = utils.dict(dict(linestyle='-', linewidth=0.8, color=axes_kw.color, alpha=0.1), **grid_kw)
        # gridminor_kw = utils.dict(dict(linestyle=':', linewidth=0.5, color=axes_kw.color, alpha=0.05), **gridminor_kw)
        # xtickminor_kw = utils.dict(xtick_kw, length=xtick_kw.length/2) # dependent on above
        # ytickminor_kw = utils.dict(ytick_kw, length=ytick_kw.length/2)
        # # ...geographic (set on declaration)
        # coastlines_kw = utils.dict(dict(linewidth=axes_kw.linewidth, color=axes_kw.color), **coastlines_kw)
        # gridlines_kw = utils.dict(dict(linestyle=':', linewidth=axes_kw.linewidth, color=axes_kw.color, alpha=0.2), **gridlines_kw)
        # # ...text (set with plt.setp, set_label_text, or on creation)
        # ticklabels_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **ticklabels_kw)
        # label_kw = utils.dict(dict(size=8, weight='normal', color=axes_kw.color), **label_kw)
        # title_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **title_kw)
        # abc_kw = utils.dict(dict(size=10, weight='normal', color=axes_kw.color), **abc_kw)
        # # ...special (set on creation)
        # legend_kw = utils.dict(dict(framealpha=0.6, fancybox=False), **legend_kw)
        # colorbar_kw = utils.dict(dict(extend='both', spacing='uniform', drawedges=False), **colorbar_kw)

        # Colorbar/legend checks
        # if (colorbar or colorbar_ax is not None) and mappable is None: # steal space from existing axis
        #     raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
        # if legend_ax is not None and objects is None:
        #     raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')
        if (self.bottomcolorbar is not None or self.rightcolorbar is not None) and mappable is None: # steal space from existing axis
            raise ValueError('To create colorbar, must supply the mappable object using mappable=<object>')
        if self.bottomlegend is not None and objects is None:
            raise ValueError('To create legend on its own axis, must supply the artists using objects=<list>')

        # Set properties
        # if ax is None:
        #     ax = self.axes # iterate through everything; apply same settings to each object
        ax = self.axes_main
        for i,a in enumerate(ax): # .flat gives row-major indexing, just like matlab for plots
            # Special axes locating
            colorbar_ax_on = False
            if a==self.bottomcolorbar: colorbar_ax_on, orientation = True, 'horizontal'
            if a==self.rightcolorbar: colorbar_ax_on, orientation = True, 'vertical'
            legend_ax_on = True if a==self.bottomlegend else False
            basemap_ax_find = np.where([a==basemap.ax for basemap in self.basemaps])[0]

            # Legend axes setup
            # ...some checks
            if legend_ax_on: # e.g. can input -1, or axes.Axes object
                legend_ax_kw = dict(legend_kw, handles=objects, loc='center')
            elif legends is not None:
                legend_ax_kw = legend_kw
            if legends is not None or legend_ax_on:
                # ...draw legend
                leg = a.legend(**legend_ax_kw)
                # ...check for certain properties only
                plt.setp(leg.get_frame(), **patch_kw)
                # ...and apply text properties, identical to ticklabel properties
                plt.setp(leg.get_texts(), **ticklabels_kw) # color should be in there as font property
                if legend_ax_on:
                    continue

            # Colorbar axes setup
            # ...some checks
            # if colorbar_ax_on:
            #     colorbar_ax_kw = dict(colorbar_kw, cax=a, use_gridspec=True, orientation=orientation) # uses space affored by this entire axis
            # elif colorbar:
            #     colorbar_ax_kw = dict(colorbar_kw, ax=a) # steals space form axis 'a', uses make_axes() to build new one
            # if colorbar or colorbar_ax_on:
            if colorbar_ax_on:
                # ...draw colorbar
                colorbar_ax_kw = dict(colorbar_kw, cax=a, use_gridspec=True, orientation=orientation) # uses space affored by this entire axis
                cb = self.colorbar(mappable, **colorbar_ax_kw)
                # ...and fix pesky white lines between levels + misalignment with border due to rasterized blocks
                cb.solids.set_rasterized(False)
                cb.solids.set_edgecolor('face')
                # ...first, make edges/dividers consistent with axis edges
                if cb.dividers is not None: plt.setp(cb.dividers, **dividers_kw)
                plt.setp(cb.outline, **patch_kw)
                # ...next, ticks formatters/locators (different from axis_setup version)
                if isinstance(clocator, mticker.Locator):
                    cb.locator = clocator
                elif clocator is not None:
                    try: iter(clocator)
                    except TypeError: clocator = utils.arange(cb.vmin, cb.vmax, clocator)
                    cb.locator = mticker.FixedLocator(clocator)
                if isinstance(cformatter, mticker.Formatter):
                    cb.formatter = cformatter
                elif cformatter is not None:
                    if isinstance(cformatter, str):
                        cb.formatter = mticker.FormatStrFormatter(cformatter) # %-style, numbers
                    else:
                        cb.formatter = mticker.FixedFormatter(cformatter)
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
                if colorbar_ax_on:
                    continue

            # Title, and abc'ing, and suptitle
            # ...title
            if titles is not None: title = titles[i]
            if title is not None: a.set_title(title, **title_kw)
            # ...suptitle
            if self.title is None and suptitle is not None:
                self.title = self.suptitle(suptitle, **title_kw)
            # ...abc
            if abc:
                a.text(0, a.title.get_position()[1], abc_convert[i].upper(),
                        transform=a.transAxes + a.titleOffsetTrans, ha='left', va='baseline',
                        **abc_kw)
                    # text method
                # plt.title(abc_convert[i].upper(), loc='left')
                #     # pyplot method
                # a.annotate(abc_convert[i].upper(), xy=(0,0), xytext=(0,1.03),
                #         xycoords='axes fraction', textcoords='axes fraction',
                #         verticalalignment='baseline', horizontalalignment='left',
                #         **abc_kw)
                #     # annotate method

            # Basemap axes setup
            if basemap_ax_find.size>0:
                # Coatlines
                m = self.basemaps[basemap_ax_find[0]]
                m.drawcoastlines(**coastlines_kw)
                # Parallels
                m.drawparallels(utils.arange(-90,90,30), **gridlines_kw)
                m.drawmeridians(utils.arange(-180,180,60), **gridlines_kw)
                continue # skip everything else
            
            # Cartopy axes setup
            # if isinstance(a, cmpl.geoaxes.GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
            #     a.add_feature(cfeature.COASTLINE, **patch_kw)
            #     plt.setp(a.outline_patch, **patch_kw)
            #     continue

            # Axes/axis properties requiring special considerations, or must be set before everything else
            # ...patch properties
            if transparent: a.patch.set_alpha(0)
            # ...spine line properties, tick line properties (should be identical between axes)
            plt.setp(list(a.spines.values()), **spines_kw)
            # ...spine positions, tick positions
            for axis, sides, spineloc, tickloc in zip((a.xaxis, a.yaxis), (('bottom','top'), ('left','right')), (xspineloc, yspineloc), (xtickloc, ytickloc)):
                # ...first, spine position
                for side in sides: 
                    if spineloc=='both':
                        a.spines[side].set_visible(True)
                    elif spineloc in sides: # make relevant spine visible
                        b = True if side==spineloc else False
                        a.spines[side].set_visible(b)
                    elif side==sides[0]: # move the left/bottom spine onto the specified location, with set_position
                        a.spines[side].set_visible(True)
                        a.spines[side].set_position(spineloc) # just gets passed to the function; options include 
                            # 'zero', 'center', and tuple with (units, location) where units can be axes, data, or outward
                    else:
                        a.spines[side].set_visible(False)
                # ...next, ticks
                if spineloc=='both' or spineloc in sides: # pass to set_ticks_position left/right/top/bottom/both
                    axis.set_ticks_position(spineloc)
                else:
                    axis.set_ticks_position(sides[0]) # will move the bottom/left spine onto zero line, etc. by default
                # ...finally, labels
                if spineloc in sides:
                    axis.set_label_position(spineloc)
                else:
                    axis.set_label_position(sides[0]) # also put label on left/bottom by default
            # ...axis scaling
            a.set_xscale(xscale)
            a.set_yscale(yscale)
            # ...axis limits, and reversal options (alternatively, supply
            # your own xlim/ylim that go from high to low)
            if xlim is None: xlim = a.get_xlim()
            if ylim is None: ylim = a.get_ylim()
            if xreverse: xlim = xlim[::-1]
            if yreverse: ylim = ylim[::-1]
            a.set_xlim(xlim)
            a.set_ylim(ylim)
            # ...and gridline stuff
            if xgrid:
                a.xaxis.grid(which='major', **grid_kw)
            if ygrid:
                a.yaxis.grid(which='major', **grid_kw)
            if xgridminor and xtickminor: # if tickminor turned off, gridminor turned on... that's dumb, ignore
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
                if isinstance(ticklocator, mticker.Locator):
                    axis.set_major_locator(ticklocator)
                elif ticklocator is not None:
                    try: iter(ticklocator)
                    except TypeError: ticklocator = utils.arange(lim[0], lim[1], ticklocator)
                    # axis.set_ticks(ticklocator)
                    axis.set_major_locator(mticker.FixedLocator(ticklocator))
                # ...next, minor tick locator
                if not tickminor:
                    axis.set_minor_locator(mticker.NullLocator()) # no minor ticks
                elif isinstance(tickminorlocator, mticker.Locator):
                    axis.set_minor_locator(tickminorlocator) # pass locator class/subclass (all documented ones subclass the parent class, mticker.Locator)
                elif tickminorlocator is not None:
                    try: iter(tickminorlocator)
                    except TypeError: tickminorlocator = utils.arange(lim[0], lim[1], tickminorlocator) # use **spacing/setp** as specified
                    axis.set_minor_locator(mticker.FixedLocator(tickminorlocator))
                # ...next, major tick formatters (enforce Null, always, for minor ticks)
                axis.set_minor_formatter(mticker.NullFormatter())
                if isinstance(tickformatter, mticker.Formatter): # use isinstance, because r"stuff" or f"stuff" **do not work**
                    axis.set_major_formatter(tickformatter)
                elif tickformatter is not None:
                    if isinstance(tickformatter, str):
                        if dates: axis.set_major_formatter(mpl.dates.DateFormatter(tickformatter)) # %-style, dates
                        else: axis.set_major_formatter(mticker.FormatStrFormatter(tickformatter)) # %-style, numbers
                    else:
                        # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
                        axis.set_major_formatter(mticker.FixedFormatter(tickformatter)) # list of strings
                plt.setp(axis.get_ticklabels(), **ticklabels_kw)
                # ...and finally, label
                if label is not None: axis.set_label_text(label)
                plt.setp(axis.get_label(), **label_kw)
                # ...and tick properties
                axis.set_tick_params(which='major', **tick_kw)
                axis.set_tick_params(which='minor', **tickminor_kw) # have length
            
        # And save, optionally (into figures folder)
        if filename is not None:
            self.save(filename, desktop=desktop)

