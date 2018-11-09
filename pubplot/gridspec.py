#!/usr/bin/env python3
import numpy as np
import matplotlib.gridspec as mgridspec
from .rcmod import rc
from .utils import dot_dict
default = lambda x,y: x if x is not None else y

# Custom settings for various journals
# Add to this throughout your career, or as standards change
# PNAS info: http://www.pnas.org/page/authors/submission
# AMS info: https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
# AGU info: https://publications.agu.org/author-resource-center/figures-faq/
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

# Function for processing input and generating necessary keyword args
def _gridspec_kwargs(array=None, rowmajor=True, # mode 1: specify everything with array
    nrows=1, ncols=1, emptycols=None, emptyrows=None, # mode 2: use convenient kwargs for simple grids
    aspect=1,    height=None, width=None,   # for controlling aspect ratio, default is control for width
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
    Parses user input and returns required ratios, spacings, and bounding
    edges for plot.
    """
    #--------------------------------------------------------------------------#
    # Parse complicated input
    #--------------------------------------------------------------------------#
    # Handle the convenience feature for specifying the desired width/spacing
    # for panels as that suitable for a colorbar or legend
    # NOTE: Ugly but this is mostly boilerplate, shouln't change much
    def _panelprops(panel, panels, colorbar, colorbars, legend, legends, width, space):
        if colorbar or colorbars:
            width = default(width, rc.subplots['cbar'])
            space = default(space, rc.subplots['labs'])
            panel, panels = colorbar, colorbars
        elif legend or legends:
            width = default(width, rc.subplots['legend'])
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
    row_offset = 0
    col_offset = 1 if leftpanels else 0

    # Apply the general defaults
    hspace = default(hspace, rc.subplots['title'])
    wspace = default(wspace, rc.subplots['inner'])
    left   = default(left,   rc.subplots['labs'])
    bottom = default(bottom, rc.subplots['labs'])
    right  = default(right,  rc.subplots['nolabs'])
    top    = default(top,    rc.subplots['title'])
    bwidth = default(bwidth, rc.subplots['cbar'])
    rwidth = default(rwidth, rc.subplots['cbar'])
    lwidth = default(lwidth, rc.subplots['cbar'])
    bspace = default(bspace, rc.subplots['labs'])
    rspace = default(rspace, rc.subplots['labs'])
    lspace = default(lspace, rc.subplots['labs'])

    # Array setup
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        order = 'C' if rowmajor else 'F' # for column major, use Fortran ordering
        array = array.reshape((nrows, ncols), order=order) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        array = array[None,:] if rowmajor else array[:,None] # interpret as single row or column
    # Empty rows/columns feature
    # TODO: Delete, outdated
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
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
    wratios = default(wratios, np.ones(ncols)/ncols) # only do this ones rows/columns figured out
    hratios = default(hratios, np.ones(nrows)/nrows)

    # Necessary arguments to reconstruct this grid
    # Can follow some of the pre-processing
    subplots_kw = dot_dict(array=array, width=width, height=height,
        hspace=hspace, wspace=wspace,
        hratios=hratios, wratios=wratios,
        nrows=nrows,   ncols=ncols,
        bottompanels=bottompanels, leftpanels=leftpanels, rightpanels=rightpanels,
        left=left,     bottom=bottom, right=right,   top=top,
        bwidth=bwidth, bspace=bspace, rwidth=rwidth, rspace=rspace, lwidth=lwidth, lspace=lspace,
        )

    # Basic figure dimension stuff
    width, height = journalsize(width, height) # if user passed width=<string>, will use that journal size
    if width is None and height is None:
        width = 5 # default behavior is use 1:1 axes, fixed width
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)

    # Prepare gridspec
    try:
        aspect = aspect[0]/aspect[1]
    except (IndexError,TypeError):
        pass # do nothing
    wspace, hspace = np.atleast_1d(wspace), np.atleast_1d(hspace)
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    wratios = np.array(wratios)/sum(wratios)
    hratios = np.array(hratios)/sum(hratios)
    aspect = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    aspect = aspect*(hratios[0]/np.mean(hratios))

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

    #--------------------------------------------------------------------------
    # Apply aspect ratio to axes, infer hspaces/wspaces/bottom/top/left/right
    #--------------------------------------------------------------------------
    # Automatic aspect ratio
    # TODO: Account for top-left axes occupying multiple subplot slots!
    bpanel_total = bwidth + bspace if bottompanels else 0
    rpanel_total = rwidth + rspace if rightpanels else 0
    lpanel_total = lwidth + lspace if leftpanels else 0
    if width is not None:
        axwidth_ave = (width - left - right - sum(wspace) - rpanel_total - lpanel_total)/ncols
    if height is not None:
        axheight_ave = (height - top - bottom - sum(hspace) - bpanel_total)/nrows
    # Fix height and top-left axes aspect ratio
    if auto_width:
        axwidth_ave = axheight_ave*aspect
        width       = axwidth_ave*ncols + left + right + sum(wspace) + rpanel_total + lpanel_total
    # Fix width and top-left axes aspect ratio
    if auto_height:
        axheight_ave = axwidth_ave/aspect
        height       = axheight_ave*nrows + top + bottom + sum(hspace) + bpanel_total
    # Figure size, and component of space belonging to main plotting area
    if axwidth_ave<0:
        raise ValueError("Not enough room for axes. Increase width, or reduce spacings 'left', 'right', or 'wspace'.")
    if axheight_ave<0:
        raise ValueError("Not enough room for axes. Increase height, or reduce spacings 'top', 'bottom', or 'hspace'.")
    axwidth_ave = (ncols*axwidth_ave + bool(rightpanels)*rwidth + bool(leftpanels)*lwidth) \
        / (ncols + bool(rightpanels) + bool(leftpanels))
    axheight_ave = (ncols*axheight_ave + bool(bottompanels)*bwidth) \
        / (nrows + bool(bottompanels))

    # Properties for outer GridSpec object
    # Make sure the 'ratios' and 'spaces' are in physical units (we cast the
    # former to physical units), easier then to add stuff as below
    wspace,  hspace  = wspace.tolist(), hspace.tolist()
    wratios, hratios = (wratios*axwidth_ave*ncols).tolist(), (hratios*axheight_ave*nrows).tolist()
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
    bottom = bottom/height
    left   = left/width
    top    = 1-top/height
    right  = 1-right/width
    wspace = [w/axwidth_ave for w in wspace]
    hspace = [h/axheight_ave for h in hspace]

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
    return figsize, array, offset, subplots_kw, gridspec_kw

# Generate custom GridSpec classes that override the GridSpecBase
# __setitem__ method and the 'base' __init__ method
def flexible_gridspec_factory(base):
    class _GridSpec(base):
        """
        Generalization of builtin matplotlib GridSpec that allows for
        subplots with *arbitrary spacing*. Accomplishes this by designating certain
        rows and columns as *empty*.

        Further accepts all spacing arguments in *inches*.

        This allows for user-specified extra spacing, and for automatic adjustment
        of spacing depending on whether labels or ticklabels are overlapping. Will
        be added to figure class as auto_adjust() method or something.
        """
        def __init__(self, nrows, ncols, **kwargs):
            # Add these as attributes; want _ratios to be self-contained, so
            # it can be invoked on already instantiated gridspec
            self._nrows_visible = nrows
            self._ncols_visible = ncols
            self._nrows = nrows*2-1
            self._ncols = ncols*2-1
            wratios, hratios, kwargs = self._ratios(**kwargs)
            return super().__init__(self._nrows, self._ncols,
                    hspace=0, wspace=0, # we implement these as invisible rows/columns
                    width_ratios=wratios,
                    height_ratios=hratios,
                    **kwargs,
                    )

        def __getitem__(self, key):
            # Interpret input
            # Note that if multiple indices are requested, e.g. gridspec[1,2], the
            # argument is passed as a tuple
            def _normalize(key, size):
                if isinstance(key, slice):
                    start, stop, _ = key.indices(size)
                    if stop > start:
                        return start, stop - 1
                else:
                    if key < 0:
                        key += size
                    if 0 <= key < size:
                        return key, key
                raise IndexError(f"Invalid index: {key} with size {size}.")
            # Get flat (row-major) numbers
            # Note SubplotSpec initialization will figure out the row/column
            # geometry of this number automatically
            # nrows, ncols = self.get_geometry()
            # nrows, ncols = self._nrows_visible, self._ncols_visible
            nrows, ncols = self._nrows, self._ncols
            nrows_visible, ncols_visible = self._nrows_visible, self._ncols_visible
            if isinstance(key, tuple):
                try:
                    k1, k2 = key
                except ValueError:
                    raise ValueError('Unrecognized subplot spec "{key}".')
                num1, num2 = np.ravel_multi_index(
                    [_normalize(k1, nrows_visible), _normalize(k2, ncols_visible)],
                    (nrows, ncols),
                    )
            else:
                num1, num2 = _normalize(key, nrows_visible * ncols_visible)
            # Adjust numbers to account for 'invisible' rows/columns
            # Think proceeding down by rows, when you move to a new column
            # that skips a 'hspace' and when you move to a new row that skips
            # a 'wspace', so just multiply by 2!
            def _adjust(n):
                if n<0:
                    n = 2*(n+1) - 1 # want -1 to stay -1, -2 becomes -3, etc.
                else:
                    n = n*2
                return n
            num1, num2 = _adjust(num1), _adjust(num2)
            return mgridspec.SubplotSpec(self, num1, num2)

        def _ratios(self, hspace=None, wspace=None, # spacing between axes
                hratios=None, wratios=None,
                height_ratios=None, width_ratios=None,
                **kwargs):
            # Allow passing x=None by user
            nrows = self._nrows_visible
            ncols = self._ncols_visible
            hratios = default(height_ratios, hratios)
            wratios = default(width_ratios, wratios)
            hspace = default(hspace, 0.05)
            wspace = default(wspace, 0.05)

            # Translate height/width ratios
            hratios = default(hratios, [1]*nrows)
            wratios = default(wratios, [1]*ncols)
            if len(hratios) != nrows:
                raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')
            if len(wratios) != ncols:
                raise ValueError(f'Got {ncols} columns, but {len(wratios)} wratios.')
            wratios = np.array(wratios)/sum(wratios)
            hratios = np.array(hratios)/sum(hratios)

            # Translate height/width spacings, implement as extra columns/rows
            try:
                if len(wspace)==1:
                    wspace = [wspace[0]]*(ncols-1) # also convert to list
            except TypeError:
                wspace = [wspace]*(ncols-1)
            try:
                if len(hspace)==1:
                    hspace = [hspace[0]]*(nrows-1)
            except TypeError:
                hspace = [hspace]*(nrows-1)
            if nrows>1 and len(hspace) != nrows-1:
                raise ValueError(f'Require {nrows-1} height spacings for {nrows} rows, got {len(hspace)}.')
            if ncols>1 and len(wspace) != ncols-1:
                raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(wspace)}.')
            wspace = np.array(wspace)*np.mean(wratios) # make relative to average axes width/height
            hspace = np.array(hspace)*np.mean(hratios)

            # Assign spacing as ratios
            wratios_final = [None]*self._ncols
            wratios_final[::2] = list(wratios)
            if self._ncols>1:
                wratios_final[1::2] = list(wspace)
            hratios_final = [None]*self._nrows
            hratios_final[::2] = list(hratios)
            if self._nrows>1:
                hratios_final[1::2] = list(hspace)
            return wratios_final, hratios_final, kwargs # bring extra kwargs back

        def update(self, **gridspec_kw):
            # Handle special hspace/wspace arguments, and just set the simple
            # left/right/top/bottom attributes
            wratios, hratios, edges_kw = self._ratios(**gridspec_kw)
            edges_kw = {key:value for key,value in edges_kw.items()
                if key not in ('nrows','ncols')} # cannot be changed
            self.set_width_ratios(wratios)
            self.set_height_ratios(hratios)
            super().update(**edges_kw) # should just be left/right/top/bottom

    return _GridSpec

# Make classes
FlexibleGridSpec = flexible_gridspec_factory(mgridspec.GridSpec)
FlexibleGridSpec.__name__ = 'FlexibleGridSpec'
FlexibleGridSpecFromSubplotSpec = flexible_gridspec_factory(mgridspec.GridSpecFromSubplotSpec)
FlexibleGridSpecFromSubplotSpec.__name__ = 'FlexibleGridSpecFromSubplotSpec'

