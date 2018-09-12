# Custom option.
# class Axes(maxes.Subplot):
#     """
#     Subclass of maxes.Subplot class, with special formatting options.
#     """
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#
#     @staticmethod
#     def from_ax(ax=None, fig=None, *args, **kwargs):
#         if ax is None:
#             if fig is None: fig = figure(facecolor='w', edgecolor='w')
#             ax = MyAxes(fig, 1, 1, 1, *args, **kwargs)
#             fig.add_axes(ax)
#         return ax

def subplots(array=None, nrows=1, ncols=1,
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True, spany=True, # custom setting; share axis labels for axes with same xmin/ymin extents
            # what about bottom vs. top xlabels, left vs. right ylabels?
        aspect=None, height=None, width=5, # for controlling aspect ratio; default is control for width
        hspace=0.5, wspace=0.5, height_ratios=None, width_ratios=None, # options for gridspec (instead of gridpsec_kw)
            # spacing betwen axes, in inches
        left=0.5, bottom=0.5, right=0.2, top=0.3,
            # edges of plot for GridSpec object, in inches
        bottomcolorbar=False, rightcolorbar=False, bottomlegend=False, # for legend/colorbar referencing multiple subplots
        bottompanel_width=0.2, rightpanel_width=0.2, # width for different panels
        colorbars=None, legends=None, # integer in array corresponding to colorbar/legend position
        package='basemap', projection=None, **projection_kw): # for projections; can be 'basemap' or 'cartopy'
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
        * Start stepping away from SubplotSpec/GridSpec, and make axes from scratch.
        ...should then be able to seamlessly make differing wspaes/hspaces.
    
    NOTE that, because I usually format all subplots with .format() method, shared axes should already
    end up with the same axis limits/scaling/majorlocators/minorlocators... the sharex and sharey detection algorithm
    really is just to get INSTRUCTIONS to make the ticklabels and axis labels INVISIBLE for certain axes.
    """
    
    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.reshape(np.arange(1,nrows*ncols+1), (nrows, ncols))
    array = np.array(array)
    for i in array.flat:
        if type(i) is not int:
            raise TypeError('Must specify subplot locations with integers.')
    nrows, ncols = array.shape[0], array.shape[1]
    
    # Projection setup stuff
    # ...for both basemap/cartopy, determine default aspect ratio
    if projection is not None and aspect is None:
        aspect = 2 # try this
        # aspects = (yrange[:,1]-yrange[:,0])/(xrange[:,1]-xrange[:,0])
        # if not np.all(aspects[0]==aspects): # should ammend to accept width/height ratios
        #     raise ValueError('If each subplot is map projection, their aspect ratios should match.')
        # aspect = .5*aspects[0] # .6 is good bet
    elif aspect is None:
        aspect = 1.5 # try this as a default; square plots (1) look too tall to me
    # ...cartopy stuff; get crs instance, and create dictionary to add to add_subplot calls
    if projection is not None and package=='cartopy':
        crs_dict = dict(
                rectilinear=ccrs.PlateCarree,
                pcarree=ccrs.PlateCarree,
                cyl=ccrs.PlateCarree, # matches basemap
                platecarree=ccrs.PlateCarree, 
                robinson=ccrs.Robinson,
                stereo=ccrs.Stereographic,
                stereographic=ccrs.Stereographic,
                moll=ccrs.Mollweide,
                mollweide=ccrs.Mollweide,
                aeqd=ccrs.AzimuthalEquidistant,
                np=ccrs.NorthPolarStereo,
                sp=ccrs.SouthPolarStereo,
                )
        cartopy_kw = dict(projection=crs_dict[projection](**projection_kw))
    else:
        cartopy_kw = dict()

    # Fixed aspect, infer required figure size [try wspace=0.2 for gs with cols=4, axes gs[0,:2] and gs[0,2:] 
    # vs. gs with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ])
    # FORMULA: aspect = ((width - left - right - (ncol-1)*wspace)/ncol) / ((height - top - bottom - (nrow-1)*hspace)/nrow)
    # ...set figure aspect ratio, fixing height
    if height is not None:
        axheight_ave_nopanel = (height - top - bottom - (nrows-1)*hspace - int(bottompanel)*(bottompanel_width + bottom))/nrows
        if axheight_ave_nopanel<0: raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axwidth_ave_nopanel = (axheight_ave_nopanel*aspect)
        width = axwidth_ave_nopanel*ncols + left + right + (ncols-1)*wspace + int(rightpanel)*(rightpanel_width + left) 
    # ...set figure aspect ratio, fixing width
    else:
        axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - int(rightpanel)*(rightpanel_width + right))/ncols
        if axwidth_ave_nopanel<0: raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axheight_ave_nopanel = (axwidth_ave_nopanel/aspect)
        height = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + int(bottompanel)*(bottompanel_width + bottom)
    # ...figure size, and some other stuff
    axwidth_total = ncols*axwidth_ave_nopanel + (ncols-1)*wspace
    axheight_total = nrows*axheight_ave_nopanel + (nrows-1)*hspace
    figsize = (width, height)

    # Convert hspace/wspace/edges to INCHES units
    width_ratios_outer = np.array([axwidth_total, rightpanel_width])/(axheight_total + rightpanel_width) if rightpanel else None
    height_ratios_outer = np.array([axheight_total, bottompanel_width])/(axheight_total + rightpanel_width) if bottompanel else None
    wspace_outer = right/((axwidth_total + rightpanel_width)/2) if rightpanel else None
    hspace_outer = bottom/((axheight_total + bottompanel_width)/2) if bottompanel else None
    wspace = wspace/axwidth_ave_nopanel
    hspace = hspace/axheight_ave_nopanel
    left, right = left/width, 1-right/width
    bottom, top = bottom/height, 1-top/height

    # Get default width/height ratios
    # ...width
    if width_ratios is None: width_ratios = np.ones(ncols)/ncols
    else: width_ratios = np.array(width_ratios)/sum(width_ratios)
    # ...height
    if height_ratios is None: height_ratios = np.ones(nrows)/nrows
    else: height_ratios = np.array(height_ratios)/sum(height_ratios)

    # Or, set up with subplotspec
    fig = plt.figure(FigureClass=Figure, figsize=figsize)
    if bottompanel or rightpanel:
        if rightpanel: right = 1-left # should be identical; planned for that so far
        GS = mgridspec.GridSpec(nrows=1+int(bottompanel), ncols=1+int(rightpanel), left=left, bottom=bottom, right=right, top=top,
                wspace=wspace_outer, hspace=hspace_outer, width_ratios=width_ratios_outer, height_ratios=height_ratios_outer) # set wspace/hspace to match the top/bottom spaces
        gs = mgridspec.GridSpecFromSubplotSpec(nrows=nrows, ncols=ncols, subplot_spec=GS[0,0],
                wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)
        # print(wspace_outer, hspace_outer, width_ratios_outer, height_ratios_outer)
    else:
        gs = mgridspec.GridSpec(nrows=nrows+int(bottompanel), ncols=ncols+int(rightpanel), # add row/column if panel exists
                left=left, bottom=bottom, right=right, top=top,
                wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)
    
    # Preliminary stuff for getting shared axes
    # # ...colorbars
    # cbar_ids = [np.where(array==-1)]
    # if colorbars is not None:
    #     try: iter(colorbars)
    #     except TypeError:
    #         colorbars = (colorbars,)
    #     # get the colorbar locations
    #     cbar_ids = [np.where(array==i) for i in colorbars]
    #     # label the colorbars
    #     for c in colorbars:
    #         array[array==colorbars] = -1
    # # ...legends
    # if legends is not None:
    #     try: iter(legends)
    #     except TypeError:
    #         legends = (legends,)
    #     # get legends locations
    #     leg_ids = [np.where(array==i) for i in legends]
    #     array[array==legends] = -2
    # ...basemaps
    # ...axes
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty, -1 for colorbar, -2 for legend
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    num_axes = len(axes_ids)
    
    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # ...get spanning x, y
    if spanx:
    def xlabelspan(self, cax, ax):
        # Much easier to make this into method, instead of having to prepare for it
        # cax is central axes (we place ylabel on that one)
        xmin = min(a.get_position().intervalx[0] for a in ax)
        xmax = max(a.get_position().intervalx[1] for a in ax)
        for a in ax:
            if a!=cax: a.xaxis.label.set_visible(False)
        cax.xaxis.label.set_transform(mtransforms.blended_transform_factory(
                self.transFigure, mtransforms.IdentityTransform(), # specify x, y transform
                ))
        cax.xaxis.label.set_position(((xmin+xmax)/2, 0))

    def ylabelspan(self, cax, ax):
        # Same as above, pretty much
        ymin = min(a.get_position().intervaly[0] for a in ax)
        ymax = max(a.get_position().intervaly[1] for a in ax)
        for a in ax:
            if a!=cax: a.yaxis.label.set_visible(False)
        cax.yaxis.label.set_transform(mtransforms.blended_transform_factory(
                mtransforms.IdentityTransform(), self.transFigure # specify x, y transform
                ))
        cax.yaxis.label.set_position((0, (ymin+ymax)/2))

    # ...get shared x
    if sharex:
        xgroups_base, xgroups, grouped = [], [], []
        for i in range(num_axes):
            matches = np.all(xrange[i:i+1,:]==xrange, axis=1) # broadcasting rules work here
            group = np.where(matches)[0]
            if i not in grouped and group.size>1:
                xgroups += [group] # add ndarray of ids to list
                xgroups_base += [group[np.argmax(yrange[group, 1])]] # bottom-most axis with shared x, for group
            grouped += [*group] # bookkeeping; record ids that have been grouped already
    # ...get shared y
    if sharey:
        ygroups_base, ygroups, grouped = [], [], []
        for i in range(num_axes):
            matches = np.all(yrange[i:i+1,:]==yrange, axis=1) # broadcasting rules work here
            group = np.where(matches)[0]
            if i not in grouped and group.size>1:
                ygroups += [group] # add ndarray of ids to list
                ygroups_base += [group[np.argmin(xrange[group, 0])]] # left-most axis with shared y, for group
            grouped += [*group] # bookkeeping; record ids that have been grouped already
    
    # Make shared axes
    ax = num_axes*[None] # list of axes
    if sharex:
        for i in xgroups_base:
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    **cartopy_kw)
    if sharey:
        for i in ygroups_base:
            if ax[i] is None:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **cartopy_kw)

    # Make dependent axes
    # ...have to manually insert the correct axis object for sharex, sharey argument on creation with add_subplot
    # print(num_axes, array, '\nxrange:', xrange, '\nyrange:', yrange, '\nxgroups, and shared:', xgroups, '\n', xgroups_base, '\nygroups, and shared:', ygroups, '\n', ygroups_base)
    for i in range(num_axes):
        if ax[i] is None: # i.e. it's not one of the axes on which others are dependent
            # ...detect if we want to share this axis with another
            sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
            if sharex:
                igroup = np.where([i in g for g in xgroups])[0] # np.where works on lists
                sharex_ax = ax[xgroups_base[igroup[0]]]
                # if igroup.size==1:
                # sharex_ax = ax[xgroups_base[igroup[0]]]
                #     if sharex_ax is None:
                #         raise ValueError('Something went wrong; shared x axes was not already drawn.')
                # elif igroup.size>1:
                #     raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)
            if sharey:
                igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
                sharey_ax = ax[ygroups_base[igroup[0]]]
                # if igroup.size==1:
                # sharey_ax = ax[ygroups_base[igroup[0]]]
                #     if sharey_ax is None:
                #         raise ValueError('Something went wrong; shared x axes was not already drawn.')
                # elif igroup.size>1:
                #     raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)
            # ...draw axis, and add to list
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], sharex=sharex_ax, sharey=sharey_ax,
                    **cartopy_kw)
            # ...hide tick labels (not default behavior for manual sharex, sharey use)
            if sharex_ax is not None:
                for t in ax[i].xaxis.get_ticklabels(): t.set_visible(False)
                ax[i].xaxis.get_label().set_visible(False)
                # plt.setp(ax[i].xaxis.get_ticklabels(), visible=False)
            if sharey_ax is not None:
                for t in ax[i].yaxis.get_ticklabels(): t.set_visible(False)
                ax[i].yaxis.get_label().set_visible(False)
                # plt.setp(ax[i].yaxis.get_ticklabels(), visible=False)
                
    # Extra stuff
    # ...create panel axes (at end)
    if bottomcolorbar:
        try: iter(bottompanel_range)
        except TypeError:
            fig.axes_colorbar = fig.add_subplot(GS[1,0])
        else:
            fig.axes_colorbar = fig.add_subplot(GS[1,0])
        fig.axes_colorbar_orientation = 'horizontal'
    if rightcolorbar:
        fig.axes_colorbar = fig.add_subplot(GS[0,1])
        fig.axes_colorbar_orientation = 'vertical'
    if bottomlegend:
        fig.axes_legend = fig.add_subplot(GS[1,0])
    # ...set basemap
    if projection is not None and package=='basemap':
        # ....make objects, and add to figure
        m = np.array(num_axes*[None])
        for i in range(num_axes):
            m[i] = mbasemap.Basemap(projection=projection, ax=ax[i], **projection_kw)
        fig.basemaps = m
        # ...if singleton, remove from list
        if len(m)==1: m = m[0]
    # ...if singleton, remove from list
    if len(ax)==1: ax = ax[0]

    # Finally, return result
    print('Axes grid constructed.')
    if projection is not None and package=='basemap':
        return fig, ax, m
    else:
        return fig, ax
