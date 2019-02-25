#!/usr/bin/env python3
"""
Module for creating special ProPlot figures with special ProPlot axes.
"""
import re
import numpy as np
# import io
# from contextlib import redirect_stdout
# Local modules, projection sand formatters and stuff
from .rcmod import rc
from .utils import _default, ic
from . import gridspec, figure, axes
import functools
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------#
# Miscellaneous helper functions
#------------------------------------------------------------------------------#
def figure(*args, **kwargs):
    """
    Alias for `subplots`.
    """
    return subplots(*args, **kwargs)

def close():
    """
    Close all figures stored in memory. Note this does not delete
    rendered figures in an iPython notebook.
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
    Magical class that iterates through items and calls respective
    method (or retrieves respective attribute) on each one. For example,
    ``axs.format(color='r')`` colors all axes spines and tick marks red.

    When using `subplots`, the list of axes returned is an instance of
    `axes_list`.
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
            @functools.wraps(values[0])
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

# Function for processing input and generating necessary keyword args
def _subplots_kwargs(nrows, ncols, rowmajor=True,
    aspect=1,    figsize=None, # for controlling aspect ratio, default is control for width
    width=None,  height=None, axwidth=None, axheight=None, journal=None,
    hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
    left=None,   bottom=None, right=None,   top=None,     # spaces around edge of main plotting area, in inches
    bwidth=None, bspace=None, rwidth=None, rspace=None, lwidth=None, lspace=None, # default to no space between panels
    bottompanel=False, bottompanels=False, # bottompanelrows=1, # optionally draw extra rows
    rightpanel=False,  rightpanels=False,  # rightpanelcols=1,
    leftpanel=False,   leftpanels=False,   # leftpanelcols=1,
    bottomcolorbar=False, bottomcolorbars=False, bottomlegend=False, bottomlegends=False, # convenient aliases that change default features
    rightcolorbar=False,  rightcolorbars=False,  rightlegend=False,  rightlegends=False,
    leftcolorbar=False,   leftcolorbars=False,   leftlegend=False,   leftlegends=False
    ):
    """
    Handle complex keyword args and aliases thereof, apply rc configuration
    defaults.
    """
    # Handle the convenience feature for specifying the desired width/spacing
    # for panels as that suitable for a colorbar or legend
    # NOTE: Ugly but this is mostly boilerplate, shouln't change much
    def _panelprops(panel, panels, colorbar, colorbars, legend, legends, width, space):
        if colorbar or colorbars:
            width = _default(width, rc['gridspec.cbar'])
            space = _default(space, rc['gridspec.xlab'])
            panel, panels = colorbar, colorbars
        elif legend or legends:
            width = _default(width, rc['gridspec.legend'])
            space = _default(space, 0)
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
    # Need to do this after number of rows/columns figured out
    wratios = np.atleast_1d(_default(wratios, 1))
    hratios = np.atleast_1d(_default(hratios, 1))
    hspace = np.atleast_1d(_default(hspace, rc['gridspec.title']))
    wspace = np.atleast_1d(_default(wspace, rc['gridspec.inner']))
    if len(wratios)==1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios)==1:
        hratios = np.repeat(hratios, (nrows,))
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    left   = units(_default(left,   rc['gridspec.ylab']))
    bottom = units(_default(bottom, rc['gridspec.xlab']))
    right  = units(_default(right,  rc['gridspec.nolab']))
    top    = units(_default(top,    rc['gridspec.title']))
    bwidth = units(_default(bwidth, rc['gridspec.cbar']))
    rwidth = units(_default(rwidth, rc['gridspec.cbar']))
    lwidth = units(_default(lwidth, rc['gridspec.cbar']))
    bspace = units(_default(bspace, rc['gridspec.xlab']))
    rspace = units(_default(rspace, rc['gridspec.ylab']))
    lspace = units(_default(lspace, rc['gridspec.ylab']))

    # Determine figure size
    if journal:
        if width or height or axwidth or axheight or figsize:
            raise ValueError('Argument conflict: Specify only a journal size, or the figure dimensions, not both.')
        width, height = journal_size(journal) # if user passed width=<string>, will use that journal size
    if not figsize:
        figsize = (width, height)
    width, height = figsize
    width  = units(width, error=False)
    height = units(height, error=False)

    # If width and height are not fixed, determine necessary width/height to
    # preserve the aspect ratio of specified plot
    auto_both = (width is None and height is None)
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    auto_neither = (width is not None and height is not None)
    bpanel_space = bwidth + bspace if bottompanels else 0
    rpanel_space = rwidth + rspace if rightpanels else 0
    lpanel_space = lwidth + lspace if leftpanels else 0
    try:
        aspect = aspect[0]/aspect[1]
    except (IndexError,TypeError):
        pass # do nothing
    aspect_fixed = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    aspect_fixed = aspect*(hratios[0]/np.mean(hratios))
    # Determine average axes widths/heights
    # Default behavior: axes average 2.0 inches wide
    if auto_width or auto_neither:
        axheight_ave = (height - top - bottom - sum(hspace) - bpanel_space)/nrows
    if auto_height or auto_neither:
        axwidth_ave = (width - left - right - sum(wspace) - rpanel_space - lpanel_space)/ncols
    if auto_both: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = 2.0
        if axheight is not None:
            height = axheight*nrows + top + bottom + sum(hspace) + bpanel_space
            auto_width = True
            axheight_ave = axheight
        if axwidth is not None:
            width = axwidth*ncols + left + right + sum(wspace) + rpanel_space + lpanel_space
            auto_height = True
            axwidth_ave = axwidth
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
        figsize = (width, height) # again
    # Fix height and top-left axes aspect ratio
    if auto_width:
        axwidth_ave = axheight_ave*aspect_fixed
        width       = axwidth_ave*ncols + left + right + sum(wspace) + rpanel_space + lpanel_space
    # Fix width and top-left axes aspect ratio
    if auto_height:
        axheight_ave = axwidth_ave/aspect_fixed
        height       = axheight_ave*nrows + top + bottom + sum(hspace) + bpanel_space
    # Check
    if axwidth_ave<0:
        raise ValueError(f"Not enough room for axes (would have width {axwidth_ave}). Increase width, or reduce spacings 'left', 'right', or 'wspace'.")
    if axheight_ave<0:
        raise ValueError(f"Not enough room for axes (would have height {axheight_ave}). Increase height, or reduce spacings 'top', 'bottom', or 'hspace'.")

    # Necessary arguments to reconstruct this grid
    # Can follow some of the pre-processing
    subplots_kw = _dot_dict(nrows=nrows, ncols=ncols,
        figsize=figsize, aspect=aspect,
        hspace=hspace,   wspace=wspace,
        hratios=hratios, wratios=wratios,
        bottompanels=bottompanels, leftpanels=leftpanels, rightpanels=rightpanels,
        left=left,     bottom=bottom, right=right,   top=top,
        bwidth=bwidth, bspace=bspace, rwidth=rwidth, rspace=rspace, lwidth=lwidth, lspace=lspace,
        )

    # Make sure the 'ratios' and 'spaces' are in physical units (we cast the
    # former to physical units), easier then to add stuff as below
    wspace = wspace.tolist()
    hspace = hspace.tolist()
    wratios = (ncols*axwidth_ave*(wratios/sum(wratios))).tolist()
    hratios = (nrows*axheight_ave*(hratios/sum(hratios))).tolist()

    # Now add the outer panel considerations (idea is we have panels whose
    # widths/heights are *in inches*, and only allow the main subplots and
    # figure widhts/heights to warp to preserve aspect ratio)
    nrows += int(bool(bottompanels))
    ncols += int(bool(rightpanels)) + int(bool(leftpanels))
    if bottompanels: # the 'bottom' space actually goes between subplots and panel
        hratios = hratios + [bwidth] # easy
        hspace  = hspace + [bottom]
        bottom  = bspace
    if leftpanels:
        wratios = [lwidth] + wratios
        wspace  = [left] + wspace
        left    = lspace
    if rightpanels:
        wratios = wratios + [rwidth]
        wspace  = wspace + [right]
        right   = rspace
    # Scale stuff that gridspec needs to be scaled
    # Scale the boundaries for gridspec
    # NOTE: We *no longer* scale wspace/hspace because we expect it to
    # be in same scale as axes ratios, much easier that way and no drawback really
    bottom = bottom/height
    left   = left/width
    top    = 1-top/height
    right  = 1-right/width

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    offset = (0, 1 if leftpanels else 0)
    figsize = (width, height)
    gridspec_kw = dict(
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
    return figsize, offset, subplots_kw, gridspec_kw

def subplots(array=None, ncols=1, nrows=1,
        order='C', # allow calling with subplots(array)
        emptycols=[], emptyrows=[], # obsolete?
        tight=None, auto_adjust=True,
        rcreset=True, silent=True, # arguments for figure instantiation
        span=None, # bulk apply to x/y axes
        share=None, # bulk apply to x/y axes
        spanx=1,  spany=1,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        sharex=3, sharey=3, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        innerpanels={}, innercolorbars={}, innerpanels_kw={}, innercolorbars_kw={},
        basemap=False, proj={}, projection={}, proj_kw={}, projection_kw={},
        **kwargs):
    # See `proplot.axes.XYAxes.format`, `proplot.axes.CartopyAxes.format`, and
    # `proplot.axes.BasemapAxes.format` for formatting your figure.
    """
    Analagous to `matplotlib.pyplot.subplots`. Create a figure with a single
    axes or arbitrary grids of axes (any of which can be map projections),
    and optional "panels" along axes or figure edges.

    Instead of separating many ProPlot features into their own functions (e.g.
    a `generate_panel` function), we specify them with keyword arguments
    right when the figure is declared.

    The reason for this approach is we want to build a
    *static figure "scaffolding"* before plotting anything. This
    way, ProPlot can tightly control axes aspect ratios and panel/colorbar
    widths. See the `auto_adjust` keyword arg for details.

    The parameters are sorted into the rough sections subplot grid
    specifications, figure and subplot sizes, axis sharing,
    inner panels, outer panels, map projections, and other settings.

    Parameters
    ----------
    ncols, nrows : int, optional
        Number of columns, rows. Ignored if `array` is not ``None``.
        Use these arguments for simpler subplot grids.
    order : {'C', 'F'}, optional
        Whether subplots are numbered in column-major (``'C'``) or row-major
        (``'F'``) order. Analagous to `numpy.array` ordering. This controls
        the order axes appear in the `axs` list, and the order of a-b-c
        labelling when using `~proplot.axes.BaseAxes.format` with ``abc=True``.
    array : None or array-like of int, optional
        2-dimensional array specifying complex grid of subplots. Think of
        this array as a "picture" of your figure. For example, the array
        ``[[1, 1], [2, 3]]`` creates one long subplot in the top row, two
        smaller subplots in the bottom row.

        Integers must range from 1 to the number of plots. For example,
        ``[[1, 4]]`` is invalid.

        ``0`` indicates an empty space. For example, ``[[1, 1, 1], [2, 0, 3]]``
        creates one long subplot in the top row with two subplots in the bottom
        row separated by a space.
    emptyrows, emptycols : list of int, optional
        Row, column numbers (starting from 1) that you want empty. Generally
        this is used with `ncols` and `nrows`; with `array`, you can
        just use zeros.

    width, height : float or str, optional
        The figure width and height. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units` -- for example,
        ``width="10cm"`` creates a 10cm wide figure.
    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    journal : None or str, optional
        Conform figure width (and height, if specified) to academic journal
        standards. See `~proplot.gridspec.journal_size` for details.
    axwidth, axheight : float or str, optional
        Sets the average width, height of your axes. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.

        These arguments are convenient where you don't care about the figure
        dimensions and just want your axes to have enough "room".
    aspect : float or length-2 list of floats, optional
        The (average) axes aspect ratio, in numeric form (width divided by
        height) or as (width, height) tuple. If you do not provide
        the `hratios` or `wratios` keyword args, all axes will have
        identical aspect ratios.
    hspace, wspace, hratios, wratios, left, right, top, bottom : optional
        Passed to `~proplot.gridspec.FlexibleGridSpecBase`.

    spanx, spany, span : bool or {1, 0}, optional
        Whether to use "spanning" axis labels for the *x* axis, *y* axis, or
        both axes. When ``True`` or ``1``, the axis label for the leftmost
        (*y*) or bottommost (*x*) subplot is **centered** on that column or
        row (i.e. it "spans" the column or row).

        This labels multiple subplot axes with one label, reducing
        redundancy of information.
    sharex, sharey, share : {3, 2, 1, 0}, optional
        The "axis sharing level" for the *x* axis, *y* axis, or both
        axes. Options are as follows:

            0. No axis sharing.
            1. Only draw *axis label* on the leftmost (*y*) or
               bottommost (*x*) column or row of subplots. Axis tick labels
               still appear on every subplot.
            2. As in 1, but forces the axis limits to be identical. Axis
               tick labels still appear on every subplot.
            3. As in 2, but only show the *axis tick labels* on the
               leftmost (*y*) or bottommost (*x*) column or row of subplots.

        This feature can considerably reduce redundancy of information
        in your figure.

    innerpanels : str or dict-like, optional
        Specify which axes should have inner "panels".

        Panels are stored on the ``leftpanel``, ``rightpanel``,
        ``bottompanel`` and ``toppanel`` attributes on the axes instance.
        You can also use the aliases ``lpanel``, ``rpanel``, ``bpanel``,
        or ``tpanel``.

        The argument is interpreted as follows:

            * If str, panels are drawn on the same side for all subplots.
              String should contain any of the characters ``'l'`` (left panel),
              ``'r'`` (right panel), ``'t'`` (top panel), or ``'b'`` (bottom panel).
              For example, ``'rt'`` indicates a right and top panel.
            * If dict-like, panels can be drawn on different sides for
              different subplots. For example, `innerpanels={1:'r', (2,3):'l'}`
              indicates that we want to draw a panel on the right side of
              subplot number 1, on the left side of subplots 2 and 3.

    innercolorbars : str or dict-like, optional
        Identical to ``innerpanels``, except the default panel size is more
        appropriate for a colorbar. The panel can then be **filled** with
        a colorbar with, for example, ``ax.leftpanel.colorbar(...)``.
    innerpanels_kw : dict-like, optional
        Keyword args passed to ``~proplot.figure.Figure.panel_factory``. Controls
        the width/height and spacing of panels.

        Can be dict of properties (applies globally), or **dict of dicts** of
        properties (applies to specific properties, as with `innerpanels`).

        For example, consider a figure with 2 columns and 1 row.
        With ``{'wwidth':1}``, all left/right panels 1 inch wide, while
        with ``{1:{'wwidth':1}, 2:{'wwidth':0.5}}``, the left subplot
        panel is 1 inch wide and the right subplot panel is 0.5 inches wide.

        See ``~proplot.figure.Figure.panel_factory`` for keyword arg options.
    innercolorbars_kw
        Alias for ``innerpanels_kw``.

    bottompanel, rightpanel, leftpanel : bool, optional
        Whether to draw an "outer" panel surrounding the entire grid of
        subplots, on the bottom, right, or left sides of the figure.
        This is great for creating global colorbars and legends.
        Default is ``False``.

        Panels are stored on the ``leftpanel``, ``rightpanel``,
        and ``bottompanel`` attributes on the figure instance.
        You can also use the aliases ``lpanel``, ``rpanel``, and ``bpanel``.
    bottompanels, rightpanels, leftpanels : bool or list of int, optional
        As with the `panel` keyword args, but this allows you to
        specify **multiple** panels along each figure edge.
        Default is ``False``.

        The arguments are interpreted as follows:

            * If bool, separate panels **for each column/row of subplots**
              are drawn.
            * If list of int, can specify panels that span **contiguous
              columns/rows**. Usage is similar to the ``array`` argument.

              For example, for a figure with 3 columns, ``bottompanels=[1, 2, 2]``
              draws a panel on the bottom of the first column and spanning the
              bottom of the right 2 columns, and ``bottompanels=[0, 2, 2]``
              only draws a panel underneath the right 2 columns (as with
              `array`, the ``0`` indicates an empty space).

    bottomcolorbar, bottomcolorbars, rightcolorbar,  rightcolorbars, leftcolorbar, leftcolorbars
        Identical to the corresponding ``panel`` and ``panels`` keyword
        args, except the default panel size is more appropriate for a colorbar.

        The panel can then be **filled** with a colorbar with, for example, 
        ``fig.leftpanel.colorbar(...)``.
    bottomlegend, bottomlegends, rightlegend,  rightlegends, leftlegend, leftlegends
        Identical to the corresponding ``panel`` and ``panels`` keyword
        args, except the default panel size is more appropriate for a legend.

        The panel can then be **filled** with a legend with, for example, 
        ``fig.leftpanel.legend(...)``.
    lwidth, rwidth, bwidth, twidth : float, optional
        Width of left, right, bottom, and top panels. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.
    lspace, rspace, bspace, tspace : float, optional
        Space between the "inner" edges of the left, right, bottom, and top
        panels and the edge of the main subplot grid. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.

    projection : str or dict-like, optional
        The map projection name.

        If string, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `innerpanels`.

        For example, consider a figure with 4 columns and 1 row.
        With ``projection={1:'mercator', (2,3):'hammer'}``,
        the leftmost subplot is a Mercator projection, the middle 2 are Hammer
        projections, and the rightmost is a normal x/y axes.
    projection_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class on instantiation.

        As with `innerpanels_kw`,
        can be dict of properties (applies globally), or **dict of dicts** of
        properties (applies to specific properties, as with `innerpanels`).
    proj, proj_kw
        Aliases for `projection`, `projection_kw`.
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.

        If boolean, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `innerpanels`.

    auto_adjust : bool, optional
        Whether to automatically adjust figure size to prevent
        `~matplotlib.artist.Artist` elements from getting cut off. Also
        standardizes the figure margin size *without* messing up axes
        aspect ratios. Default is ``True``.
    tight : bool, optional
        Alias for ``auto_adjust``.
    rcreset : bool, optional
        Whether to reset all `plot.rc` settings to their default values
        once the figure is drawn.

    silent : bool, optional
        Whether to suppress print statements. Default is ``True``.

    Returns
    -------
    f : `~proplot.figure.Figure`
        The figure instance.
    axs : `~proplot.figure.axes_list`
        A special list of axes instances. See `~proplot.figure.axes_list`.

    Other parameters
    ----------------
    **kwargs
        Passed to ``~proplot.gridspec.FlexibleGridSpec``.

    Notes
    -----
    Matplotlib `~matplotlib.axes.Axes.set_aspect` option seems to behave
    strangely for some plots; for this reason we override the ``fix_aspect``
    keyword arg provided by `~mpl_toolkits.basemap.Basemap` and just draw
    the figure with appropriate aspect ratio to begin with. Otherwise, we
    get weird differently-shaped subplots that seem to make no sense.

    Shared axes will generally end up with the same axis
    limits/scaling/majorlocators/minorlocators; the ``sharex`` and ``sharey``
    detection algorithm really is just to get instructions to make the
    ticklabels/axis labels **invisible** for certain axes.

    Todo
    ----
    * Add options for e.g. ``bpanel`` keyword args.
    * Fix axes aspect ratio stuff when width/height ratios are not one!
    * Generalize axes sharing for right y-axes and top x-axes. Enable a secondary
      axes sharing mode where we *disable ticklabels and labels*, but *do not
      use the builtin sharex/sharey API*, suitable for complex map projections.
    * For spanning axes labels, right now only detect **x labels on bottom**
      and **ylabels on top**. Generalize for all subplot edges.
    """
    # Check
    sharex = _default(share, sharex)
    sharey = _default(share, sharey)
    spanx  = _default(span, spanx)
    spany  = _default(span, spany)
    if int(sharex) not in range(4) or int(sharey) not in range(4):
        raise ValueError('Axis sharing options sharex/sharey can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels).')
    # Helper functions
    translate = lambda p: {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}.get(p, p)
    auto_adjust = _default(tight, auto_adjust)
    def axes_dict(value, kw=False):
        # First build up dictionary
        # Accepts:
        # 1) 'string' or {1:'string1', (2,3):'string2'}
        if not kw:
            if not isinstance(value, dict):
                value = {range(1,num_axes+1): value}
        # 2) {'prop':value} or {1:{'prop':value1}, (2,3):{'prop':value2}}
        else:
            nested = [isinstance(value,dict) for value in value.values()]
            if not any(nested): # any([]) == False
                value = {range(1,num_axes+1): value.copy()}
            elif not all(nested):
                raise ValueError('Wut.')
        # Then unfurl wherever keys contain multiple axes numbers
        kw_out = {}
        for nums,item in value.items():
            nums = np.atleast_1d(nums)
            for num in nums.flat:
                if not kw:
                    kw_out[num-1] = item
                else:
                    kw_out[num-1] = item.copy()
        # Verify numbers
        if {*range(num_axes)} != {*kw_out.keys()}:
            raise ValueError(f'Have {num_axes} axes, but {value} has properties for axes {", ".join(str(i+1) for i in sorted(kw_out.keys()))}.')
        return kw_out

    # Array setup
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        if order not in ('C','F'): # better error message
            raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
        array = array.reshape((nrows, ncols), order=order) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        if order not in ('C','F'): # better error message
            raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
        array = array[None,:] if order=='C' else array[:,None] # interpret as single row or column
    # Empty rows/columns feature
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    if emptycols:
        emptycols = np.atleast_1d(emptycols)
        for col in emptycols.flat:
            array[:,col-1] = 0
    if emptyrows:
        emptyrows = np.atleast_1d(emptyrows)
        for row in emptyrows.flat:
            array[row-1,:] = 0
    # Enforce rule
    nums = np.unique(array[array!=0])
    num_axes = len(nums)
    if tuple(nums.flat) != tuple(range(1,num_axes+1)):
        raise ValueError('Axes numbers must span integers 1 to num_axes (i.e. cannot skip over numbers).')
    nrows = array.shape[0]
    ncols = array.shape[1]

    # Get basemap.Basemap or cartopy.CRS instances for map, and override aspec tratio
    # NOTE: Previously went to some pains (mainly for basemap, something in the
    # initialization deals with this) to only draw one projection. This is hard
    # to generalize when want different projections/kwargs, so abandon
    basemap = axes_dict(basemap, False) # package used for projection
    proj = axes_dict(projection or proj or 'xy', False) # name of projection; by default use base.XYAxes
    proj_kw = axes_dict(projection_kw or proj_kw, True) # stores cartopy/basemap arguments
    axes_kw = {num:{} for num in range(num_axes)} # stores add_subplot arguments
    for num,name in proj.items():
        # Builtin matplotlib polar axes, just use my overridden version
        if name=='polar':
            axes_kw[num]['projection'] = 'newpolar'
            if num==1:
                kwargs.update(aspect=1)
        # The default, my XYAxes projection
        elif name=='xy':
            axes_kw[num]['projection'] = 'xy'
        # Custom Basemap and Cartopy axes
        elif name:
            package = 'basemap' if basemap[num] else 'cartopy'
            instance, aspect = axes.map_projection_factory(package, name, **proj_kw[num])
            axes_kw[num].update({'projection':package, 'map_projection':instance})
            if not silent:
                print(f'Forcing aspect ratio: {aspect:.3g}')
            if num==0:
                kwargs.update(aspect=aspect)
        else:
            raise ValueError('All projection names should be declared. Wut.')

    # Create dictionary of panel toggles and settings
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    # NOTE: Internally we convert array references to 0-base here
    # Add kwargs and the 'which' arguments
    # Optionally change the default panel widths for 'colorbar' panels
    if not isinstance(innercolorbars, (dict, str)):
        raise ValueError('Must pass string of panel sides or dictionary mapping axes numbers to sides.')
    if not isinstance(innerpanels, (dict, str)):
        raise ValueError('Must pass string of panel sides or dictionary mapping axes numbers to sides.')
    innerpanels = axes_dict(innerpanels or '', False)
    innercolorbars = axes_dict(innercolorbars or '', False)
    if innercolorbars_kw:
        innerpanels_kw = innercolorbars_kw
    innerpanels_kw = axes_dict(innerpanels_kw, True)
    for kw in innerpanels_kw.values():
        kw['whichpanels'] = ''
    for num,which in innerpanels.items():
        innerpanels_kw[num]['whichpanels'] += translate(which)
    for num,which in innercolorbars.items():
        which = translate(which)
        if which:
            innerpanels_kw[num]['whichpanels'] += which
            if re.search('[bt]', which):
                kwargs['hspace'] = _default(kwargs.get('hspace',None), rc['gridspec.xlab'])
                innerpanels_kw[num]['sharex_panels'] = False
                innerpanels_kw[num]['hwidth'] = _default(innerpanels_kw[num].get('hwidth', None), rc['gridspec.cbar'])
                innerpanels_kw[num]['hspace'] = _default(innerpanels_kw[num].get('hspace', None), rc['gridspec.xlab'])
            if re.search('[lr]', which):
                kwargs['wspace'] = _default(kwargs.get('wspace',None), rc['gridspec.ylab'])
                innerpanels_kw[num]['sharey_panels'] = False
                innerpanels_kw[num]['wwidth'] = _default(innerpanels_kw[num].get('wwidth', None), rc['gridspec.cbar'])
                if 'l' in which and 'r' in which:
                    default = (rc['gridspec.ylab'], rc['gridspec.nolab'])
                elif 'l' in which:
                    default = rc['gridspec.ylab']
                else:
                    default = rc['gridspec.nolab']
                innerpanels_kw[num]['wspace'] = _default(innerpanels_kw[num].get('wspace', None), default)

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    figsize, offset, subplots_kw, gridspec_kw = \
            _subplots_kwargs(nrows, ncols, **kwargs)
    row_offset, col_offset = offset
    gs = gridspec.FlexibleGridSpec(**gridspec_kw)
    fig = plt.figure(figsize=figsize, auto_adjust=auto_adjust, rcreset=rcreset,
        gridspec=gs, subplots_kw=subplots_kw,
        FigureClass=figure.Figure,
        )

    #--------------------------------------------------------------------------
    # Manage shared axes/axes with spanning labels
    #--------------------------------------------------------------------------
    # Get some axes properties
    # Note that these locations should be **sorted** by axes id
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty
    yrange = row_offset + np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = col_offset + np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    # asdfas
    # xmin   = np.array([xy[0].min() for xy in axes_ids]) # unused
    # ymax   = np.array([xy[1].max() for xy in axes_ids])

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
        ax_kw = axes_kw[i]
        if axs[i] is not None: # already created
            continue
        if innerpanels_kw[i]['whichpanels']: # non-empty
            axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    spanx=spanx, spany=spany,
                    number=i+1, **ax_kw, **innerpanels_kw[i]) # main axes handle
        else:
            axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    spanx=spanx, spany=spany,
                    number=i+1, **ax_kw) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another. If so, get that
        # axes. Also do some error checking
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        ax_kw = axes_kw[i]
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
                axs[i]._sharex_setup(sharex_ax, sharex)
            if sharey_ax is not None and axs[i] is not sharey_ax:
                axs[i]._sharey_setup(sharey_ax, sharey)
        else:
            # Virgin axes; these are not an x base or a y base
            if innerpanels_kw[i]['whichpanels']: # non-empty
                axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1, spanx=spanx, spany=spany,
                        sharex_level=sharex, sharey_level=sharey,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw, **innerpanels_kw[i])
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1, spanx=spanx, spany=spany,
                        sharex_level=sharex, sharey_level=sharey,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw) # main axes can be a cartopy projection

    # Check that axes don't belong to multiple groups
    # This should be impossible unless my code is completely wrong...
    for ax in axs:
        for name,groups in zip(('sharex', 'sharey'), (xgroups, ygroups)):
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
            axp = fig.add_subplot(subspec, panel_side=side, invisible=True, projection='panel')
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

