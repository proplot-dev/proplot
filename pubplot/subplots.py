#!/usr/bin/env python3
import re
import numpy as np
# import io
# from contextlib import redirect_stdout
import matplotlib.pyplot as plt
# Local modules, projection sand formatters and stuff
from .rcmod import rc
from .gridspec import _gridspec_kwargs, FlexibleGridSpec
from . import utils
from . import base
from functools import wraps
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

    def __getattr__(self, attr):
        # Stealthily return dummy function that actually iterates
        # through each attribute here
        values = [getattr(ax, attr, None) for ax in self]
        if None in values:
            raise AttributeError(f"'{type(self[0])}' object has no method '{attr}'.")
        elif all(callable(value) for value in values):
            @wraps(values[0])
            def iterator(*args, **kwargs):
                ret = []
                for ax in self:
                    res = getattr(ax, attr)(*args, **kwargs)
                    if res is not None:
                        ret += [res]
                return None if not ret else ret[0] if len(ret)==1 else ret
            return iterator
        elif all(not callable(value) for value in values):
            return values[0] if len(values)==1 else values # just return the attribute list
        else:
            raise AttributeError('Mixed methods found.')

def subplots(array=None, # allow calling with subplots(array)
        tight=True, rcreset=True, silent=True, # arguments for figure instantiation
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True,  spany=True,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        ihspace=None, iwspace=None, ihwidth=None, iwwidth=None,
        innerpanels=None, innercolorbars=None,
        whichpanel=None, whichpanels='r', # same as below; list of numbers where we want subplotspecs
        maps=None, # set maps to True for maps everywhere, or to list of numbers
        basemap=False, proj=None, projection=None, proj_kw={}, projection_kw={},
        **kwargs): # for projections; can be 'basemap' or 'cartopy'
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
            ihwidth = default(ihwidth, rc.subplots['cbar'])
            ihspace  = default(ihspace, rc.subplots['labs'])
        elif re.search('[lr]', whichpanels):
            iwwidth = default(iwwidth, rc.subplots['cbar'])
            iwspace  = default(iwspace,  rc.subplots['labs'])
    # Translate whichpanels
    whichpanels = whichpanel or whichpanels # user can specify either
    translate = {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}
    whichpanels = translate.get(whichpanels, whichpanels)

    # Get basemap.Basemap or cartopy.CRS instances for map, and override aspec tratio
    map_kw = {}
    projection = proj or projection
    projection_kw = proj_kw or projection_kw
    if maps:
        wratios = np.atleast_1d(kwargs.get('wratios',None) or 1)
        hratios = np.atleast_1d(kwargs.get('hratios',None) or 1)
        if {*wratios.flat} != {1} or {*hratios.flat} != {1}:
            raise NotImplementedError('Not yet possible.')
        if projection=='polar':
            map_kw = {'projection':'newpolar'}
            kwargs.update(aspect=1)
        else:
            package = 'basemap' if basemap else 'cartopy'
            map_projection, aspect = base.map_projection_factory(package, projection, **projection_kw)
            map_kw = {'projection':package, 'map_projection':map_projection}
            if not silent:
                print(f'Forcing aspect ratio: {aspect:.3g}')
            kwargs.update(aspect=aspect)

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    figsize, array, offset, subplots_kw, gridspec_kw = _gridspec_kwargs(array=array, **kwargs)
    row_offset, col_offset = offset
    gs = FlexibleGridSpec(**gridspec_kw)
    fig = plt.figure(figsize=figsize, tight=tight, rcreset=rcreset,
        gridspec=gs, subplots_kw=subplots_kw,
        FigureClass=base.Figure,
        )

    # Find axes that use map projections
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
    panel_kw = {'whichpanels':whichpanels,
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
    yrange = row_offset + np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = col_offset + np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
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
        ax_kw = map_kw if i in maps_ids else {'projection':'xy'}
        if axs[i] is not None: # already created
            continue
        # print('base', i)
        if i in innerpanels_ids:
            axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    number=i+1, **ax_kw, **panel_kw) # main axes handle
        else:
            axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    number=i+1, **ax_kw) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another. If so, get that
        # axes. Also do some error checking
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        ax_kw = map_kw if i in maps_ids else {'projection':'xy'}
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
                        number=i+1,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw, **panel_kw)
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw) # main axes can be a cartopy projection

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

    #--------------------------------------------------------------------------#
    # Create panel axes
    #--------------------------------------------------------------------------#
    def _paneladd(name, panels):
        if not panels:
            return
        axsp = []
        side = re.sub('^(.*)panel$', r'\1', name)
        for n in np.unique(panels).flat:
            offset = row_offset if side in ('left','right') else col_offset
            idx, = np.where(panels==n)
            idx = slice(offset + min(idx), offset + max(idx) + 1)
            if side=='right':
                subspec = gs[idx,-1]
            elif side=='left':
                subspec = gs[idx,0]
            elif side=='bottom':
                subspec = gs[-1,idx]
            axp = fig.add_subplot(subspec, panelside=side, invisible=True, projection='panel')
            axsp += [axp]
        setattr(fig, name, axes_list(axsp))
    _paneladd('bottompanel', subplots_kw.bottompanels)
    _paneladd('rightpanel',  subplots_kw.rightpanels)
    _paneladd('leftpanel',   subplots_kw.leftpanels)

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

