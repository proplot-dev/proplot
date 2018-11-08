#!/usr/bin/env python3
import re
import numpy as np
# import io
# from contextlib import redirect_stdout
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
# Local modules, projection sand formatters and stuff
from .gridspec import FlexibleGridSpec, FlexibleGridSpecFromSubplotSpec
from . import utils
from . import base
# Conversions
cm2in = 0.3937
mm2in = cm2in/10.0
in2cm = 1.0/cm2in
in2mm = 1.0/mm2in
default = lambda x,y: x if x is not None else y

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
class axes_list(list):
    """
    Special list of axes, iterates through each axes and calls respective
    method on each one. Will return a list of each return value.
    """
    def __repr__(self):
        # Make clear that this is no ordinary list
        return 'axes_list(' + super().__repr__() + ')'

    def __getitem__(self, key):
        # Return an axes_list version of the slice, or just the axes
        axs = list.__getitem__(self, key)
        if isinstance(key,slice): # i.e. returns a list
            axs = axes_list(axs)
        return axs

    def __getattribute__(self, attr):
        # Stealthily return dummy function that actually iterates
        # through each attribute here
        for ax in self:
            if not hasattr(ax,attr) or not callable(getattr(ax,attr)):
                raise AttributeError(f"'{type(ax)}' object has no method '{attr}'.")
        def iterator(*args, **kwargs):
            ret = []
            for ax in self:
                ret += [getattr(ax,attr)(*args, **kwargs)]
            return ret
        return iterator

def subplots(array=None, rowmajor=True, # mode 1: specify everything with array
        nrows=1, ncols=1, emptycols=None, emptyrows=None, # mode 2: use convenient kwargs for simple grids
        tight=False,  # whether to set up tight bbox from gridspec object
        silent=True, # how much stuff to print
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True,  spany=True,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        aspect=1,    height=None, width=None,   # for controlling aspect ratio, default is control for width
        hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
        left=None,   bottom=None, right=None,   top=None,     # spaces around edge of main plotting area, in inches
        bwidth=None, bspace=None, rwidth=None, rspace=None, lwidth=None, lspace=None, # default to no space between panels
        bottompanel=False, bottompanels=False, bottompanelrows=1, # optionally draw extra rows
        rightpanel=False,  rightpanels=False,  # rightpanelcols=1,
        leftpanel=False,   leftpanels=False,   # leftpanelcols=1,
        bottomcolorbar=False, bottomcolorbars=False, bottomlegend=False, bottomlegends=False, # convenient aliases that change default features
        rightcolorbar=False,  rightcolorbars=False,  rightlegend=False,  rightlegends=False,
        leftcolorbar=False,   leftcolorbars=False,   leftlegend=False,   leftlegends=False,
        space_title = 0.4,   # extra space for title/suptitle
        space_inner = 0.2,   # just have ticks, no labeels
        space_legend = 0.25, # default legend space (bottom of figure)
        space_cbar = 0.17,   # default colorbar width
        space_labs = 0.5,    # default space wherever we expect tick and axis labels (a bit large if axis has no negative numbers/minus sign tick labels)
        space_nolabs = 0.15, # only ticks
        innerpanels=None, innercolorbars=None, whichpanel=None, whichpanels='r', # same as below; list of numbers where we want subplotspecs
        ihspace=None, iwspace=None, ihwidth=None, iwwidth=None,
        maps=None, # set maps to True for maps everywhere, or to list of numbers
        basemap=False, projection=None, projection_dict={}, **projection_kwargs): # for projections; can be 'basemap' or 'cartopy'
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
    # Handle the convenience feature for specifying width of inner panels
    # to be that suitable for a colorbar
    # Argument can be innercolorbars='b' or innercolorbars=True (means all)
    # TODO: Add this, but not that important
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
    # Translate whichpanels
    whichpanels = whichpanel or whichpanels # user can specify either
    translate = {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}
    whichpanels = translate.get(whichpanels, whichpanels)

    # Handle the convenience feature for specifying the desired width/spacing
    # for panels as that suitable for a colorbar or legend
    # NOTE: Ugly but this is mostly boilerplate, shouln't change much
    def _panelprops(panel, panels, colorbar, colorbars, legend, legends, width, space):
        if colorbar or colorbars:
            width = default(width, space_cbar)
            space = default(space, space_labs)
            panel, panels = colorbar, colorbars
        elif legend or legends:
            width = default(width, space_legend)
            space = default(space, 0)
            panel, panels = legend, legends
        return panel, panels, width, space
    rightpanel, rightpanels, rwidth, rspace, = _panelprops(
        rightpanel, rightpanels, rightcolorbar, rightcolorbars,
        rightlegend, rightlegends, rwidth, rspace)
    leftpanel, leftpanels, lwidth, lspace = _panelprops(
        leftpanel, leftpanels, leftcolorbar, leftcolorbars,
        leftlegend, leftlegends, lwidth, lspace)
    bottompanel, bottompanels, bwidth, bspace = _panelprops(
        bottompanel, bottompanels, bottomcolorbar, bottomcolorbars,
        bottomlegend, bottomlegends, bwidth, bspace)

    # Handle the convenience feature for generating one panel per row/column
    # and one single panel for all rows/columns
    def _parse(panel, panels, nmax):
        if panel: # one spanning panel
            panels = [1]*nmax
        elif panels not in (None,False): # can't test truthiness, want user to be allowed to pass numpy vector!
            try:
                panels = list(panels)
            except TypeError:
                panels = [*range(nmax)] # pass True to make panel for each column
        return panels
    bottompanels = _parse(bottompanel, bottompanels, ncols)
    rightpanels  = _parse(rightpanel,  rightpanels,  nrows)
    leftpanels   = _parse(leftpanel,   leftpanels,   nrows)

    # Apply the general defaults
    hspace = default(hspace, space_title)
    wspace = default(wspace, space_inner)
    left   = default(left, space_labs)
    bottom = default(bottom, space_labs)
    right  = default(right, space_nolabs)
    top    = default(top, space_title)
    bwidth = default(bwidth, space_cbar)
    rwidth = default(rwidth, space_cbar)
    lwidth = default(lwidth, space_cbar)
    bspace = default(bspace, space_labs)
    rspace = default(rspace, space_labs)
    lspace = default(lspace, space_labs)
    # Basic figure dimension stuff
    width, height = journalsize(width, height) # if user passed width=<string>, will use that journal size
    if width is None and height is None:
        width = 5 # default behavior is use 1:1 axes, fixed width
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    try:
        aspect = aspect[0]/aspect[1]
    except TypeError:
        pass # do nothing
    # Prepare gridspec
    main_rows = [0, nrows-1]
    main_cols = [1, ncols] if leftpanels else [0, ncols-1]
    if wratios is not None:
        wratios = np.array(wratios)/sum(wratios)
    else:
        wratios = np.ones(ncols)/ncols
    if hratios is not None:
        hratios = np.array(hratios)/sum(hratios)
    else:
        hratios = np.ones(nrows)/nrows
    aspect = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    aspect = aspect*(hratios[0]/np.mean(hratios))
    # Save some properties so they can be added to the figure
    gridprops = dict(
        left=left,     bottom=bottom, right=right,   top=top,
        bwidth=bwidth, bspace=bspace, rwidth=rwidth, rspace=rspace,
        auto_width=auto_width, auto_height=auto_height,
        )


    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        order = 'C' if rowmajor else 'F' # for column major, use Fortran ordering
        array = array.reshape((nrows, ncols), order=order) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        array = array[None,:] if rowmajor else array[:,None] # interpret as single row or column
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    # # Enforce consistent numbering; row-major increasing from 1 every time
    # # a new axes is encountered; e.g. [[1,2],[1,3],[1,4],[1,4]]
    # number = 1
    # newarray = np.zeros(array.shape)
    # for row in newarray.shape[0]:
    #     for col in newarray.shape[1]:
    #         if array[row,col] not in array.flat: # not already declared
    #             newarray[array==array[row,col]] = number
    #             number += 1
    # array = newarray

    # Include functionality to declare certian rows/columns empty, if user requested it
    # Useful if user wants to include extra space between some columns greater than the
    # spaces between axes in other columns
    if emptycols is not None:
        emptycols = np.atleast_1d(emptycols)
        for col in emptycols.flat:
            array[:,col-1] = 0
    if emptyrows is not None:
        emptyrows = np.atleast_1d(emptyrows)
        for row in emptyrows.flat:
            array[row-1,:] = 0
    nrows = array.shape[0]
    ncols = array.shape[1]

    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # Preliminary stuff
    # Note that these locations should be **sorted** by axes id
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty
    yrange = main_rows[0] + np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = main_cols[0] + np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    xmin   = np.array([xy[0].min() for xy in axes_ids])
    ymax   = np.array([xy[1].max() for xy in axes_ids])
    num_axes = len(axes_ids)

    # Get basemap.Basemap or cartopy.CRS instances for map
    # projection axes
    map_kwargs = {}
    if maps and ({*wratios.flat} != {1} or {*hratios.flat} != {1}):
        raise NotImplementedError('Not yet possible.')
    if not maps and projection_kwargs:
        raise ValueError(f'Unknown kwargs: {", ".join(projection_kwargs.keys())}. If you want to create maps, you must change the "maps" kwarg.')
    if maps:
        if projection=='polar':
            map_kwargs = {'projection':'newpolar'}
        else:
            proj_kwargs = {**projection_kwargs, **projection_dict}
            package = 'basemap' if basemap else 'cartopy'
            map_projection, aspect = base.map_projection_factory(package, projection, **proj_kwargs)
            map_kwargs = {'projection':package, 'map_projection':map_projection}

    #--------------------------------------------------------------------------
    # Apply aspect ratio to axes, infer hspaces/wspaces/bottom/top/left/right
    #--------------------------------------------------------------------------
    # Automatic aspect ratio
    # TODO: Account for top-left axes occupying multiple subplot slots!
    # * Try wspace=0.2 for gs with cols=4, axes gs[0,:2] and gs[0,2:] vs. gs
    #   with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ
    # * Formula: aspect = ((width - left - right - (ncol-1)*wspace)/ncol)
    #                   / ((height - top - bottom - (nrow-1)*hspace)/nrow)
    bpanel_total = bwidth + bspace if bottompanels else 0
    rpanel_total = rwidth + rspace if rightpanels else 0
    lpanel_total = lwidth + lspace if leftpanels else 0
    if width is not None:
        axwidth_ave = (width - left - right - (ncols-1)*np.mean(wspace) - rpanel_total - lpanel_total)/ncols
    if height is not None:
        axheight_ave = (height - top - bottom - (nrows-1)*np.mean(hspace) - bpanel_total)/nrows
    # Fix height and top-left axes aspect ratio
    if auto_width:
        axwidth_ave = axheight_ave*aspect
        width       = axwidth_ave*ncols + left + right + (ncols-1)*np.mean(wspace) + rpanel_total + lpanel_total
    # Fix width and top-left axes aspect ratio
    if auto_height:
        axheight_ave = axwidth_ave/aspect
        height       = axheight_ave*nrows + top + bottom + (nrows-1)*np.mean(hspace) + bpanel_total
    # Figure size, and component of space belonging to main plotting area
    if axwidth_ave<0:
        raise ValueError("Not enough room for axes. Increase width, or reduce spacings 'left', 'right', or 'wspace'.")
    if axheight_ave<0:
        raise ValueError("Not enough room for axes. Increase height, or reduce spacings 'top', 'bottom', or 'hspace'.")
    axwidth_total_nopanel  = ncols*axwidth_ave + (ncols-1)*np.mean(wspace)
    axheight_total_nopanel = nrows*axheight_ave + (nrows-1)*np.mean(hspace)

    # Properties for inner GridSpec object
    # NOTE: Actually do *not* have to make wspace/hspace vectors to pass to
    # our FlexibleGridSpec, but needed for further processing the *panel axes*
    # down the line below! Remove this later!
    wspace, hspace = np.atleast_1d(wspace), np.atleast_1d(hspace)
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    wspace_orig = wspace
    hspace_orig = hspace
    wspace = np.atleast_1d(wspace)/axwidth_ave
    hspace = np.atleast_1d(hspace)/axheight_ave

    # Properties for outer GridSpec object, need borders in fractional units
    # Ratios and stuff
    nrows += int(bool(bottompanels))
    ncols += int(bool(rightpanels)) + int(bool(leftpanels))
    wspace, hspace = wspace.tolist(), hspace.tolist()
    wratios, hratios = (wratios*axwidth_total_nopanel).tolist(), (hratios*axheight_total_nopanel).tolist()
    if bottompanels: # the 'bottom' space actually goes between subplots and panel
        hratios = hratios + [bwidth] # easy
        hspace = hspace + [bottom/((bwidth + axheight_total_nopanel)/2)]
        bottom = bspace
    if leftpanels:
        wratios = [lwidth] + wratios
        wspace = [left/((lwidth + axwidth_total_nopanel)/2)] + wspace
        left = lspace
    if rightpanels:
        wratios = wratios + [rwidth]
        rspace = rspace + [right/((rwidth + axwidth_total_nopanel)/2)]
        right = rspace
    bottom = bottom/height
    left   = left/width
    top    = 1-top/height
    right  = 1-right/width

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    # GS = mgridspec.GridSpec(
    GS = FlexibleGridSpec(
            nrows         = nrows,
            ncols         = ncols,
            left          = left,
            bottom        = bottom,
            right         = right, # unique spacing considerations
            top           = top, # so far no panels allowed here
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = wratios,
            height_ratios = hratios,
            ) # set wspace/hspace to match the top/bottom spaces
    fig = plt.figure(gridspec=GS, gridprops=gridprops, figsize=(width, height),
        FigureClass=base.Figure,
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
            axs[i] = fig.panel_factory(GS[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                    number=i+1, **ax_kwargs, **panel_kwargs) # main axes handle
        else:
            axs[i] = fig.add_subplot(GS[slice(*yrange[i,:]), slice(*xrange[i,:])],
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
                axs[i] = fig.panel_factory(GS[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        width, height, number=i+1,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kwargs, **panel_kwargs)
            else:
                axs[i] = fig.add_subplot(GS[slice(*yrange[i,:]), slice(*xrange[i,:])],
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
    # This draws
    for ax in axs:
        ax._panel_setup()

    # Check that axes don't belong to multiple groups
    # This should be impossible unless my code is completely wrong...
    for ax in axs:
        for name,groups in zip(('sharex', 'sharey', 'spanx', 'spany'), (xgroups, ygroups, xgroups_span, ygroups_span)):
            if sum(ax in group for group in xgroups)>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple {name} groups.')

    # Create panel axes
    def _paneladd(name, panels):
        if not panels:
            return
        axsp = []
        side = re.sub('^(.*)panel$', r'\1', name)
        for n in np.unique(panels).flat:
            offset = main_rows[0] if side in ('left','right') else main_cols[0]
            idx, = np.where(panels==n)
            idx = slice(offset + min(idx), offset + max(idx) + 1)
            if side=='right': # NOTE: Indexing -1 causes weird error!
                subspec = GS[idx,main_cols[1]+1]
            elif side=='left':
                subspec = GS[idx,0]
            elif side=='bottom':
                subspec = GS[main_rows[1]+1,idx]
            axp = fig.add_subplot(subspec, panelside=side, invisible=True, projection='panel')
            axsp += [axp]
        setattr(fig, name, axes_list(axsp))
    _paneladd('bottompanel', bottompanels)
    _paneladd('rightpanel',  rightpanels)
    _paneladd('leftpanel',   leftpanels)

    #--------------------------------------------------------------------------
    # Return results
    # Will square singleton arrays
    #--------------------------------------------------------------------------
    if not silent:
        print('Figure setup complete.')
    # if len(axs)==1:
    #     axs = axs[0]
    # return fig, axs
    return fig, axes_list(axs)

