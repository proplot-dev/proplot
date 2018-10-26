#!/usr/bin/env python3
import re
import numpy as np
# import io
# from contextlib import redirect_stdout
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
# Local modules, projection sand formatters and stuff
from . import utils
from . import base

#------------------------------------------------------------------------------#
# Miscellaneous helper functions
#------------------------------------------------------------------------------#
def figure(*args, **kwargs):
    """
    Simple alias for 'subplots', perhaps more intuitive.
    """
    return subplots(*args, **kwargs)
def close():
    """
    Close all figures 'open' in memory. This does not delete images printed
    in an ipython notebook; those are rendered versions of the abstract figure objects.
    """
    plt.close('all') # easy peasy
def show():
    """
    Show all figures.
    """
    plt.show()

#------------------------------------------------------------------------------#
# Custom settings for various journals
# Add to this throughout your career, or as standards change
# PNAS info: http://www.pnas.org/page/authors/submission
# AMS info: https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
# AGU info: https://publications.agu.org/author-resource-center/figures-faq/
#------------------------------------------------------------------------------#
def journalsize(width, height):
    # User wants to define their own size
    if type(width) is not str:
        return width, height
    # Determine automatically
    cm2in = 0.3937
    mm2in = 0.03937
    table = {
        'pnas1': 8.7*cm2in,
        'pnas2': 11.4*cm2in,
        'pnas3': 17.8*cm2in,
        'ams1': 3.2,
        'ams2': 4.5,
        'ams3': 5.5,
        'ams4': 6.5,
        'agu1': (95*mm2in, 115*mm2in),
        'agu2': (190*mm2in, 115*mm2in),
        'agu3': (95*mm2in, 230*mm2in),
        'agu4': (190*mm2in, 230*mm2in),
        }
    value = table.get(width, None)
    if value is None:
        raise ValueError(f'Unknown journal figure width specifier "{width}". ' +
                          'Options are: ' + ', '.join(table.keys()))
    # Output width, and optionally also specify height
    if utils.isnumber(value):
        width = value
    else:
        width, height = value
    return width, height

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
def share():
    """
    Checks group of shared axes, applies shared axes settings to the ones
    that are already drawn.
    """
    pass

def span():
    """
    Checks group of spanning axes, applies shared axes settings to the ones
    that are already drawn.
    """
    pass

def subplots(array=None, rowmajor=True, # mode 1: specify everything with array
        nrows=1, ncols=1, emptycols=None, emptyrows=None, # mode 2: use convenient kwargs for simple grids
        tight=False,  # whether to set up tight bbox from gridspec object
        silent=True, # how much stuff to print
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True,  spany=True,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        aspect=1,    height=None, width=None,   # for controlling aspect ratio, default is control for width
        hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
        left=None,   bottom=None, right=None,   top=None,     # spaces around edge of main plotting area, in inches
        bwidth=None, bspace=None, rwidth=None, rspace=None, # default to no space between panels
        bottompanel=False,    bottompanels=False,    rightpanel=False,    rightpanels=False,    bottompanelrows=1,  rightpanelcols=1,    # optionally draw extra rows
        bottomcolorbar=False, bottomcolorbars=False, rightcolorbar=False, rightcolorbars=False, bottomlegend=False, bottomlegends=False,
        space_title = 0.4,   # extra space for title/suptitle
        space_inner = 0.2,   # just have ticks, no labeels
        space_legend = 0.25, # default legend space (bottom of figure)
        space_cbar = 0.17,   # default colorbar width
        space_labs = 0.5,    # default space wherever we expect tick and axis labels (a bit large if axis has no negative numbers/minus sign tick labels)
        space_nolabs = 0.15, # only ticks
        innerpanels=None, innercolorbars=None, whichpanel=None, whichpanels='r', # same as below; list of numbers where we want subplotspecs
        ihspace=None, iwspace=None, ihwidth=None, iwwidth=None,
        maps=None, # set maps to True for maps everywhere, or to list of numbers
        package='basemap', projection=None, projection_dict={}, **projection_kwargs): # for projections; can be 'basemap' or 'cartopy'
    """
    Summary
    -------
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes objects.

    Details
    -------
    * Easiest way to create subplots is with nrows=1 and ncols=1. If you want extra space
      between a row or column, specify the row/column number that you want to be 'empty' with
      emptyrows=row/emptycolumn=column, and adjust wratios/hratios for the desired width of that space.
    * For more complicated plots, can pass e.g. array=[[1,2,3,4],[0,5,5,0]] to create a grid
      of 4 plots on the top, single plot spanning the middle 2-columns on the bottom, and empty
      spaces where the 0 appears.
    * Use bottompanel/bottompanels to make several or multiple panels on the bottom
      that can be populated with multiple colorbars/legend; bottompanels=True will
      just make one 'space' for every column, and bottompanels=[1,1,2] for example will
      make a panel spanning the first two columns, then a single panel for the final column.
      This will add a bottompanel attribute to the figure; can index that attribute if there
      are multiple places for colorbars/legend.
    * Use rightpanel/rightpanels in the same way.
    * Create extra panels *within* a grid of subplots (e.g. a 2x2 grid, each of which has a teeny
      panel attached) innerpanels=True, and optionally filter to subplot numbers with whichpanels=[list];
      then use ihspace/iwspace/ihwidth/iwwidth to control the separation and widths of the subpanels.
    * Initialize cartopy plots with package='basemap' or package='cartopy'. Can control which plots
      we want to be maps with maps=True (everything) or maps=[numbers] (the specified subplot numbers).

    Notes
    -----
    * Matplotlib set_aspect option seems to behave strangely on some plots (trend-plots from
        SST paper); for this reason we override the fix_aspect option provided by basemap and
        just draw figure with appropriate aspect ratio to begin with. Otherwise get weird
        differently-shaped subplots that seem to make no sense.
    * Shared axes will generally end up with the same axis limits/scaling/majorlocators/minorlocators;
        the sharex and sharey detection algorithm really is just to get instructions to make the
        ticklabels/axis labels invisible for certain axes.

    Todo
    ----
    * For spanning axes labels, right now only detect **x labels on bottom**
        and **ylabels on top**; generalize for all subplot edges.
    * Figure size should be constrained by the dimensions of the axes, not vice
        versa; might make things easier.
    """
    # Translate whichpanels
    whichpanels = whichpanel or whichpanels # user can specify either
    translate = {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}
    whichpanels = translate.get(whichpanels, whichpanels)
    # First translate arguments and set context dependent defaults
    # This is convenience feature to get more sensible default spacing
    default = lambda x,y: x if x is not None else y
    if rightcolorbar or rightcolorbars:
        rwidth = default(rwidth, space_cbar)
        rspace = default(rspace, space_labs)
        rightpanel  = rightcolorbar
        rightpanels = rightcolorbars
    if bottomcolorbar or bottomcolorbars:
        bwidth = default(bwidth, space_cbar)
        bspace = default(bspace, space_labs)
        bottompanel  = bottomcolorbar
        bottompanels = bottomcolorbars
    if bottomlegend or bottomlegends:
        bwidth = default(bwidth, space_legend)
        bspace = default(bspace, 0)
        bottompanel  = bottomlegend
        bottompanels = bottomlegends
    if innercolorbars is not None:
        innerpanels = innercolorbars
        if re.search('[bt]', whichpanels):
            ihwidth = default(ihwidth, space_cbar)
            hspace  = default(hspace, space_labs)
            if 'b' in whichpanels: bottom = default(bottom, space_labs)
            if 't' in whichpanels: top = default(top, space_labs)
        elif re.search('[lr]', whichpanels):
            iwwidth = default(iwwidth, space_cbar)
            wspace  = default(wspace,  space_labs)
            if 'l' in whichpanels: left = default(left, space_labs)
            if 'r' in whichpanels: right = default(right, space_labs)
    # Next the general defaults
    hspace = default(hspace, space_title)
    wspace = default(wspace, space_inner)
    left   = default(left, space_labs)
    bottom = default(bottom, space_labs)
    right  = default(right, space_nolabs)
    top    = default(top, space_title)
    bwidth = default(bwidth, space_cbar)
    rwidth = default(rwidth, space_cbar)
    bspace = default(bspace, space_labs)
    rspace = default(rspace, space_labs)
    # spaces around edge of main plotting area, in inches
    # Apply remaining defaults
    # Fill in some basic required settings
    if package not in ['basemap','cartopy']:
        raise ValueError("Plotting package must be one of basemap or cartopy.")
    width, height = journalsize(width, height) # if user passed width=<string>, will use that journal size
    if width is None and height is None: # at least one must have value
        width = 5
    if width is not None and height is not None: # specify exact size
        aspect = width/height
    if wratios is not None:
        aspect = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    if hratios is not None:
        aspect = aspect*(hratios[0]/np.mean(hratios))
    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        array = array.reshape((nrows, ncols)) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        array = array[None,:] if rowmajor else array[:,None] # interpret as single row or column
    elif not rowmajor:
        array = np.reshape(array.flatten(), (array.shape[0],array.shape[1]), order='F') # make column major
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    # Include functionality to declare certian rows/columns empty, if user requested it
    # Useful if user wants to include extra space between some columns greater than the
    # spaces between axes in other columns
    if emptycols is not None:
        if utils.isnumber(emptycols):
            emptycols = [emptycols]
        for col in emptycols:
            array[:,col-1] = 0
    if emptyrows is not None:
        if utils.isnumber(emptyrows):
            emptyrows = [emptyrows]
        for row in emptyrows:
            array[row-1,:] = 0
    nrows = array.shape[0]
    ncols = array.shape[1]
    # # Enforce consistent numbering; row-major increasing from 1 every time
    # # a new axes is encountered; e.g. [[1,2],[1,3],[1,4],[1,4]]
    # # Maybe ignore for now, but this makes sure the output axes list is
    # # always in same row-major order, not random user-input order
    # number = 1
    # newarray = np.zeros(array.shape)
    # for row in newarray.shape[0]:
    #     for col in newarray.shape[1]:
    #         if array[row,col] not in array.flat: # not already declared
    #             newarray[array==array[row,col]] = number
    #             number += 1

    #--------------------------------------------------------------------------
    # Get basemap.Basemap or cartopy.CRS instances for map
    # projection axes
    #--------------------------------------------------------------------------
    map_kwargs = {}
    if maps and (wratios is not None or hratios is not None):
        raise NotImplementedError('Not yet possible.')
    if not maps and projection_kwargs:
        raise ValueError(f'Unknown kwargs: {", ".join(projection_kwargs.keys())}. If you want to create maps, you must change the "maps" kwarg.')
    if maps:
        if projection=='polar':
            map_kwargs = {'projection':'newpolar'}
        else:
            proj_kwargs = {**projection_kwargs, **projection_dict}
            map_projection, aspect = base.map_projection_factory(package, projection, **proj_kwargs)
            map_kwargs = {'projection':package, 'map_projection':map_projection}

    #--------------------------------------------------------------------------
    # DelayedPanel considerations; keep things general for future imporvements
    #--------------------------------------------------------------------------
    # Spacings due to panel considerations
    if bottompanel or bottompanels:
        bottom_extra_axes = bwidth
        bottom_extra_hspace = bspace
    else:
        bottom_extra_axes, bottom_extra_hspace = 0, 0
    bottom_extra = bottom_extra_axes + bottom_extra_hspace

    # And right spacings
    if rightpanel or rightpanels:
        right_extra_axes   = rwidth
        right_extra_wspace = rspace
    else:
        right_extra_axes, right_extra_wspace = 0, 0
    right_extra = right_extra_axes + right_extra_wspace

    #--------------------------------------------------------------------------
    # Apply aspect ratio to axes, infer hspaces/wspaces/bottom/top/left/right
    #--------------------------------------------------------------------------
    # Fixed aspect, infer required figure size [try wspace=0.2 for gs with cols=4, axes gs[0,:2] and gs[0,2:] 
    # vs. gs with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ])
    # FORMULA: aspect = ((width - left - right - (ncol-1)*wspace)/ncol) / ((height - top - bottom - (nrow-1)*hspace)/nrow)
    # Fixing height, determine figure aspect ratio
    if height is not None:
        axheight_ave_nopanel = (height - top - bottom - (nrows-1)*hspace - bottom_extra)/nrows
        if axheight_ave_nopanel<0:
            raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axwidth_ave_nopanel = axheight_ave_nopanel*aspect
        width               = axwidth_ave_nopanel*ncols + left + right + (ncols-1)*wspace + right_extra

    # Fixing width, determine aspect ratio
    else:
        # print(width, left, right, ncols, wspace, right_extra, ncols)
        # print(width, left, right, ncols, wspace, right_extra)
        axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - right_extra)/ncols
        if axwidth_ave_nopanel<0:
            raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axheight_ave_nopanel = axwidth_ave_nopanel/aspect
        height               = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + bottom_extra

    # Figure size, and some other stuff
    # The "total" axes width include hspace/wspace 
    axwidth_total  = ncols*axwidth_ave_nopanel + (ncols-1)*wspace
    axheight_total = nrows*axheight_ave_nopanel + (nrows-1)*hspace
    figsize        = (width, height)

    # Properties for outer GridSpec object, need borders in fractional units
    ileft, ibottom, iright, itop = left, bottom, right, top
    left = left/width
    top  = 1-top/height
    if bottom_extra_axes:
        hspace_outer        = bottom/((axheight_total + bottom_extra_axes)/2) # same for main axes and bottom panel, but 'bottom'
        height_ratios_outer = np.array([axheight_total, bottom_extra_axes])/(axheight_total + bottom_extra_axes)
        bottom              = bottom_extra_hspace/height
    else:
        hspace_outer        = 0
        height_ratios_outer = np.array([1])
        bottom              = bottom/height
    if right_extra_axes:
        width_ratios_outer = np.array([axwidth_total, right_extra_axes])/(axwidth_total + right_extra_axes)
        wspace_outer       = right/((axwidth_total + right_extra_axes)/2) # want space between main axes and panel to be 'right'
        right              = 1-right_extra_wspace/width
    else:
        width_ratios_outer = np.array([1])
        wspace_outer       = 0
        right              = 1-right/width
    # Properties for inner GridSpec object
    wspace_orig = wspace
    hspace_orig = hspace
    wspace = wspace/axwidth_ave_nopanel
    hspace = hspace/axheight_ave_nopanel
    if wratios is not None:
        wratios = np.array(wratios)/sum(wratios)
    else:
        wratios = np.ones(ncols)/ncols
    if hratios is not None:
        hratios = np.array(hratios)/sum(hratios)
    else:
        hratios = np.ones(nrows)/nrows

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    GS = mgridspec.GridSpec(
            nrows         = 1+int(bottom_extra_axes>0),
            ncols         = 1+int(right_extra_axes>0),
            left          = left,
            bottom        = bottom,
            right         = right, # unique spacing considerations
            top           = top, # so far no panels allowed here
            wspace        = wspace_outer,
            hspace        = hspace_outer,
            width_ratios  = width_ratios_outer,
            height_ratios = height_ratios_outer,
            ) # set wspace/hspace to match the top/bottom spaces
    # Create axes within the 'main' plotting area
    # Will draw individual axes using this GridSpec object later
    gs = mgridspec.GridSpecFromSubplotSpec(
            nrows         = nrows,
            ncols         = ncols,
            subplot_spec  = GS[0,0],
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = wratios,
            height_ratios = hratios,
            )
    # Create figure
    fig = plt.figure(figsize=figsize, FigureClass=base.Figure,
        left=ileft,    bottom=ibottom, right=iright,  top=itop,
        bwidth=bwidth, bspace=bspace,  rwidth=rwidth, rspace=rspace,
        height=height, width=width,    gridspec=GS
        )

    #--------------------------------------------------------------------------#
    # Selective attributes; make some axes map axes, and others simple
    # panels and stuff
    #--------------------------------------------------------------------------#
    # Find axes that use map projections
    # TODO: Implement this, not currently used! Implementation will be similar to
    # innerpanels_ids below.
    if maps is not None:
        if maps is True:
            maps = [*np.unique(array)]
        if utils.isnumber(maps):
            maps = [maps] # just want a single map
        else:
            maps = [*maps] # force into list, not array
        maps_ids = [i for i,n in enumerate(np.unique(array[array>0])) if n in maps]
    else:
        maps_ids = []
    # Find axes that have inner panels
    panel_kwargs = {'whichpanels':whichpanels,
        'hspace':ihspace, 'wspace':iwspace,
        'hwidth':ihwidth, 'wwidth':iwwidth,
        }
    if innerpanels is not None:
        if innerpanels is True:
            innerpanels = [*np.unique(array)]
        if utils.isnumber(innerpanels):
            innerpanels = [innerpanels]
        else:
            innerpanels = [*innerpanels] # force into list, not array
        innerpanels_ids = [i for i,n in enumerate(np.unique(array[array>0])) if n in innerpanels]
    else:
        innerpanels_ids = []

    #--------------------------------------------------------------------------
    # Manage shared axes/axes with spanning labels
    #--------------------------------------------------------------------------
    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # Preliminary stuff
    # Note that these locations should be **sorted** by axes id
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    xmin   = np.array([xy[0].min() for xy in axes_ids])
    ymax   = np.array([xy[1].max() for xy in axes_ids])
    num_axes = len(axes_ids)

    # Find pairs with edges on same gridspec
    xgroups_span_base, xgroups_span, grouped = [], [], []
    if spanx:
        for i in range(num_axes):
            matching_axes = np.where(xmin[i]==xmin)[0]
            if i not in grouped and matching_axes.size>1:
                # Add the axes that is farthest down to control the x-label
                xgroups_span      += [matching_axes] # add ndarray of ids to list
                xgroups_span_base += [matching_axes[np.argmin(xmin[matching_axes])]]
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    ygroups_span_base, ygroups_span, grouped = [], [], []
    if spany:
        for i in range(num_axes):
            matching_axes = np.where(ymax[i]==ymax)[0]
            if i not in grouped and matching_axes.size>1:
                # Add the axes that is farthest left
                ygroups_span      += [matching_axes] # add ndarray of ids to list
                ygroups_span_base += [matching_axes[np.argmax(ymax[matching_axes])]]
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already

    # Shared axes: generate list of base axes-dependent axes pairs
    # That is, find where the minimum-maximum gridspec extent in 'x' for a
    # given axes matches the minimum-maximum gridspec extent for a base axes
    xgroups_base, xgroups_sorted, xgroups, grouped = [], [], [], []
    if sharex:
        for i in range(num_axes): # axes now have pseudo-numbers from 0 to num_axes-1
            matches       = (xrange[i,:]==xrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0] # gives ID number of matching_axes, from 0 to num_axes-1
            if i not in grouped and matching_axes.size>1:
                # Find all axes that have the same gridspec 'x' extents
                xgroups      += [matching_axes]
                # Get bottom-most axis with shared x; should be single number
                # xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
                xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
                # Sorted group
                xgroups_sorted += [matching_axes[np.argsort(yrange[matching_axes,1])[::-1]]] # bottom-most axes is first
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    ygroups_base, ygroups_sorted, ygroups, grouped = [], [], [], []
    if sharey:
        for i in range(num_axes):
            matches       = (yrange[i,:]==yrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups      += [matching_axes]
                ygroups_base += [matching_axes[np.argmin(xrange[matching_axes,0])]] # left-most axis with shared y, for matching_axes
                ygroups_sorted += [matching_axes[np.argsort(xrange[matching_axes,0])]] # left-most axis is first
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already

    #--------------------------------------------------------------------------
    # Draw axes
    # TODO: Need to configure to automatically determine 'base' axes based on
    # what has already been drawn. Not critical but would be nice.
    # TODO: Need to do something similar for the spanning axes. Also will
    # allow label to be set on any of the axes, but when this happens, will
    # set the label on the 'base' spanning axes.
    #--------------------------------------------------------------------------
    # Base axes; to be shared with other axes as ._sharex, ._sharey attributes
    axs = num_axes*[None] # list of axes
    allgroups_base = []
    if sharex:
        allgroups_base += xgroups_base
    if sharey:
        allgroups_base += ygroups_base
    for i in allgroups_base:
        ax_kwargs = map_kwargs if i in maps_ids else {'projection':'xy'}
        if axs[i] is not None: # already created
            continue
        # print('base', i)
        if i in innerpanels_ids:
            axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                    number=i+1, **ax_kwargs, **panel_kwargs) # main axes handle
        else:
            axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    number=i+1, **ax_kwargs) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another. If so, get that
        # axes. Also do some error checking
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        ax_kwargs = map_kwargs if i in maps_ids else {'projection':'xy'}
        if sharex:
            igroup = np.where([i in g for g in xgroups])[0]
            if igroup.size==1:
                sharex_ax = axs[xgroups_base[igroup[0]]]
                if sharex_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
        if sharey:
            igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
            if igroup.size==1:
                sharey_ax = axs[ygroups_base[igroup[0]]]
                if sharey_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')

        # Draw axes, and add to list
        if axs[i] is not None:
            # Axes is a *base* and has already been drawn, but might still
            # have shared axes (e.g. is bottom-axes of three-column plot
            # and we want it to share the leftmost y-axis)
            if sharex_ax is not None and axs[i] is not sharex_ax:
                axs[i]._sharex_setup(sharex_ax)
            if sharey_ax is not None and axs[i] is not sharey_ax:
                axs[i]._sharey_setup(sharey_ax)
        else:
            # Virgin axes; these are not an x base or a y base
            if i in innerpanels_ids:
                axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        width, height, number=i+1,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kwargs, **panel_kwargs)
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kwargs) # main axes can be a cartopy projection

    # Spanning axes; allow xlabels/ylabels to span them
    # TODO: The panel setup step has to be called *after* spanx/spany/sharex/sharey
    # setup is called, but we also need to call spanx/spany setup only *after*
    # every axes is drawn. Instead, then, we call the panel setup method below
    if spanx and len(xgroups_span)>0:
        for g,b in zip(xgroups_span, xgroups_span_base):
            axs[b]._spanx_setup([axs[i] for i in g])

    if spany and len(ygroups_span)>0:
        for g,b in zip(ygroups_span, ygroups_span_base):
            axs[b]._spany_setup([axs[i] for i in g])

    # Call panel setup after everything is done
    for ax in axs:
        ax._panel_setup()

    # Check that axes don't belong to multiple groups
    # This should be impossible unless my code is completely wrong...
    for ax in axs:
        for name,groups in zip(('sharex', 'sharey', 'spanx', 'spany'), (xgroups, ygroups, xgroups_span, ygroups_span)):
            if sum(ax in group for group in xgroups)>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple {name} groups.')

    #---------------------------------------------------------------------------
    # Create panel axes
    # TODO: Fix lformat and cformat methods to apply to SubplotSpec objects
    # so they can be drawn immediately; can pass the figure handle needed to draw
    # the axes with add_subplot using a defaulted kwarg fig=fig
    #---------------------------------------------------------------------------
    # First the bottompanel options
    if bottompanel: # one spanning panel
        bottompanels = [1]*ncols
    elif bottompanels: # more flexible spec
        if utils.isnumber(bottompanels):
            bottompanels = [*range(ncols)] # pass True to make panel for each column
    if bottompanels: # non-empty
        # First get the number of columns requested and their width 
        # ratios to align with main plotting area
        num = bottompanels[0]
        npanel = len(bottompanels)
        widths_orig = [axwidth_ave_nopanel*ncols*wratio for wratio in wratios] # assumes wratios normalized
        widths_new = [widths_orig[0]] # start with this
        for p,panel in enumerate(bottompanels):
            if p==0:
                continue
            newnum = panel
            if newnum==num: # same as previous number
                widths_new[-1] += widths_orig[p] + wspace_orig
                npanel -= 1 # we want one less panel object
            else:
                widths_new += [widths_orig[p]] # add to width
            num = newnum
        rratios = [width/sum(widths_new) for width in widths_new] # just normalize it
        rspace = wspace_orig/(sum(widths_new)/npanel) # divide by new average width
        # Now create new GridSpec and add each of its
        # SubplotSpecs to the figure instance
        P = mgridspec.GridSpecFromSubplotSpec(
                nrows         = bottompanelrows,
                ncols         = npanel,
                subplot_spec  = GS[1,0],
                wspace        = rspace, # same as above
                width_ratios  = rratios,
                )
        panels = [fig.add_subplot(P[0,i], side='bottom', erase=True, projection='panel') for i in range(npanel)]
        if bottompanel:
            panels = panels[0] # no indexing if user specified single panel, but does mean indexing if specified bottompanels=True with single column grid
        fig.bottompanel = panels

    # Next the rightpanel options (very similar)
    if rightpanel:
        rightpanels = [1]*nrows
    elif rightpanels:
        if utils.isnumber(rightpanels):
            rightpanels = [*range(nrows)] # pass True to make panel for each row
    if rightpanels: # non-empty
        # First get the number of columns requested and their width 
        # ratios to align with main plotting area
        num = rightpanels[0]
        npanel = len(rightpanels)
        heights_orig = [axheight_ave_nopanel*nrows*hratio for hratio in hratios] # assumes wratios normalized
        heights_new = [heights_orig[0]] # start with this
        for p,panel in enumerate(rightpanels):
            if p==0:
                continue
            newnum = panel
            if newnum==num: # same as previous number
                heights_new[-1] += heights_orig[p] + hspace_orig
                npanel -= 1 # we want one less panel object
            else:
                heights_new += [heights_orig[p]] # add to width
            num = newnum
        rratios = [height/sum(heights_new) for height in heights_new] # just normalize it
        rspace = hspace_orig/(sum(heights_new)/npanel) # divide by new average width
        # Now create new GridSpec and add each of its
        # SubplotSpecs to the figure instance
        P = mgridspec.GridSpecFromSubplotSpec(
                ncols         = rightpanelcols,
                nrows         = npanel,
                subplot_spec  = GS[0,1],
                hspace        = rspace,
                height_ratios = rratios,
                )
        panels = [fig.add_subplot(P[0,i], side='right', erase=True, projection='panel') for i in range(npanel)]
        if rightpanel:
            panels = panels[0]
        fig.rightpanel = panels

    #--------------------------------------------------------------------------
    # Return results
    # Will square singleton arrays
    #--------------------------------------------------------------------------
    if not silent:
        print('Figure setup complete.')
    if len(axs)==1:
        axs = axs[0]
    return fig, axs

