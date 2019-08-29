#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures and axes.
The `subplots` function is all you'll need to directly use here.
It returns a `Figure` instance and an `axes_grid` container of
`~proplot.axes.Axes` axes, whose positions are controlled by the
`FlexibleGridSpec` class.

.. raw:: html

   <h1>Developer notes</h1>

Matplotlib's `~matplotlib.pyplot.subplots` returns a numpy `~numpy.ndarray`
of axes. ProPlot's `subplots` returns an `axes_grid` of axes, which is a `list`
subclass invented to implement some new features. But why does `axes_grid`
subclass `list` instead of `~numpy.ndarray`? Two reasons.

1. ProPlot's `subplots` function is meant to work with *arbitrary* arangements of
   subplots that may span *multiple* rows and columns -- i.e. the subplots don't
   necessarily fit cleanly as entries in a 2D matrix. `axes_grid` is really a 1D
   list, but permits 2D indexing *just in case* the user *happened* to draw a
   clean 2D matrix of subplots. 1D indexing of `axes_grid` (e.g.  ``axs[0]``)
   always returns a single axes, while 1D indexing of a 2D `~numpy.ndarray` of
   axes would just return a row of axes.

2. Further, since `~axes_grid.__getattr__` iterates through each axes in the
   container and returns a *single item* if the container is *singleton*, it
   doesn't really matter whether the `~proplot.subplots.subplots` return value is
   an axes or a singleton `axes_grid` of axes. This allows ProPlot to *always*
   return an `axes_grid` when `~proplot.subplots.subplots` is called. Matplotlib's
   `~matplotlib.pyplot.subplots`, by contrast, returns a single axes when one
   subplot is drawn, a 1D `~numpy.ndarray` when a single row or column is drawn,
   and a 2D `~numpy.ndarray` when a multi-row, multi-column figure is drawn.
"""
# NOTE: Importing backend causes issues with sphinx, and anyway not sure it's
# always included, so make it optional
import os
import re
import numpy as np
import functools
import warnings
import matplotlib.pyplot as plt
import matplotlib.figure as mfigure
import matplotlib.transforms as mtransforms
import matplotlib.gridspec as mgridspec
try:
    import matplotlib.backends.backend_macosx as mbackend
except ImportError:
    mbackend = None
from .rctools import rc
from .utils import _notNone, _counter, units
from . import projs, axes
__all__ = [
    'axes_grid', 'close', 'show', 'subplots', 'Figure',
    'FlexibleGridSpec', 'FlexibleGridSpecBase', 'FlexibleGridSpecFromSubplotSpec',
    ]

# Translation
_side_translate = {
    'l':'left',
    'r':'right',
    'b':'bottom',
    't':'top',
    }

#-----------------------------------------------------------------------------#
# Miscellaneous stuff
#-----------------------------------------------------------------------------#
# Wrapper functions, so user doesn't have to import pyplot
def close():
    """Alias for ``matplotlib.pyplot.close('all')``, included so you don't have
    to import `~matplotlib.pyplot`. Closes all figures stored
    in memory."""
    plt.close('all') # easy peasy

def show():
    """Alias for ``matplotlib.pyplot.show()``, included so you don't have
    to import `~matplotlib.pyplot`. Note this command should be
    unnecessary if you are doing inline iPython notebook plotting and ran the
    `~proplot.notebook.nbsetup` command."""
    plt.show()

# Helper classes
class axes_grid(list):
    """List subclass and pseudo-2D array that is used as a container for the
    list of axes returned by `subplots`, lists of figure panels, and lists of
    stacked axes panels. The shape of the array is stored in the ``shape``
    attribute. See the `~axes_grid.__getattr__` and `~axes_grid.__getitem__`
    methods for details."""
    def __init__(self, list_, n=1, order='C'):
        # Add special attributes that support 2D grids of axes
        # NOTE: The input list is always a vector *already unfurled* in row-major
        # or column-major order, and 'n' is the fastest-moving dimension size, i.e.
        # ncols for order == 'C' and nrows for order == 'F'.
        self._n = n # ncols or nrows
        self._order = order # order
        super().__init__(list_)
        self.shape = (len(self)//n, n)[::(1 if order == 'C' else -1)]

    def __repr__(self):
        return 'axes_grid(' + super().__repr__() + ')'

    def __setitem__(self, key, value):
        """Pseudo immutability, raises error."""
        raise LookupError('axes_grid is immutable.')

    def __getitem__(self, key):
        """If an integer is passed, the item is returned, and if a slice is passed,
        an `axes_grid` of the items is returned. You can also use 2D indexing,
        and the corresponding axes in the axes grid will be chosen.

        Example
        -------

        >>> import proplot as plot
        ... f, axs = plot.subplots(nrows=3, ncols=3, colorbars='b', bstack=2)
        ... axs[0] # the subplot in the top-right corner
        ... axs[3] # the first subplot in the second row
        ... axs[1,2] # the subplot in the second row, third from the left
        ... axs[:,0] # the subplots in the first column

        """
        # Allow 2D specification
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]
        if not isinstance(key, tuple): # do not expand single slice to list of integers or we get recursion! len() operator uses __getitem__!
            axlist = isinstance(key, slice)
            objs = list.__getitem__(self, key)
        elif len(key) == 2:
            axlist = any(isinstance(ikey, slice) for ikey in key)
            # Expand keys
            keys = []
            order = self._order
            for i,ikey in enumerate(key):
                if (i == 1 and order == 'C') or (i == 0 and order != 'C'):
                    n = self._n
                else:
                    n = len(self)//self._n
                if isinstance(ikey, slice):
                    start, stop, step = ikey.start, ikey.stop, ikey.step
                    if start is None:
                        start = 0
                    elif start < 0:
                        start = n + start
                    if stop is None:
                        stop = n
                    elif stop < 0:
                        stop = n + stop
                    if step is None:
                        step = 1
                    ikeys = [*range(start, stop, step)]
                else:
                    if ikey < 0:
                        ikey = n + ikey
                    ikeys = [ikey]
                keys.append(ikeys)
            # Get index pairs and get objects
            # Note that in double for loop, right loop varies fastest, so
            # e.g. axs[:,:] delvers (0,0), (0,1), ..., (0,N), (1,0), ...
            # Remember for order == 'F', axes_grid was sent a list unfurled in
            # column-major order, so we replicate row-major indexing syntax by
            # reversing the order of the keys.
            objs = []
            if self._order == 'C':
                idxs = [key0*self._n + key1 for key0 in keys[0] for key1 in keys[1]]
            else:
                idxs = [key1*self._n + key0 for key1 in keys[1] for key0 in keys[0]]
            for idx in idxs:
                objs.append(list.__getitem__(self, idx))
            if not axlist: # objs will always be length 1
                objs = objs[0]
        else:
            raise IndexError
        # Return
        if axlist:
            return axes_grid(objs)
        else:
            return objs

    def __getattr__(self, attr):
        """If the attribute is *callable*, returns a dummy function that loops
        through each identically named method, calls them in succession, and
        returns a tuple of the results. This lets you call arbitrary methods
        on multiple axes at once! If the `axes_grid` has length ``1``,
        just returns the single result. If the attribute is *not callable*,
        returns an `axes_grid` of identically named attributes for every object
        in the list.

        Example
        -------

        >>> import proplot as plot
        ... f, axs = plot.subplots(nrows=2, ncols=2, axcolorbars='b')
        ... axs.format(xticks=5) # calls "format" on all subplots in the list
        ... axs.bpanel.colorbar(m) # calls "colorbar" on all panels in the axes_grid returned by "axs.bpanel"

        """
        attrs = (*(getattr(ax, attr) for ax in self),) # may raise error
        # Empty
        if not attrs:
            def null_iterator(*args, **kwargs):
                return None
            return null_iterator
        # Panels
        if all(isinstance(_, (axes_grid, axes.Axes, axes.EmptyPanel))
            for _ in attrs):
            return axes_grid(attrs)
        # Objects
        elif not any(callable(_) for _ in attrs):
            if len(self) == 1:
                return attrs[0]
            else:
                return attrs
        # Methods
        elif all(callable(_) for _ in attrs):
            @functools.wraps(attrs[0])
            def axes_grid_iterator(*args, **kwargs):
                ret = []
                for func in attrs:
                    ret.append(func(*args, **kwargs))
                ret = (*ret,)
                if len(self) == 1:
                    return ret[0]
                elif all(res is None for res in ret):
                    return None
                else:
                    return ret
            return axes_grid_iterator
        # Mixed
        raise AttributeError(f'Found mixed types for attribute {attr!r}.')

#-----------------------------------------------------------------------------#
# Gridspec classes
#-----------------------------------------------------------------------------#
def _positem(n):
    """Account for negative indices."""
    if n < 0:
        return 2*(n+1) - 1 # want -1 to stay -1, -2 becomes -3, etc.
    else:
        return n*2

def _normalize(key, size):
    """Transform gridspec index into standardized form."""
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

class FlexibleGridSpecBase(object):
    """
    Generalization of builtin `~matplotlib.gridspec.GridSpec` that allows for
    grids with **arbitrary spacing** between rows and columns of axes.

    Accomplishes this by actually drawing ``nrows*2 + 1`` and ``ncols*2 + 1``
    `~matplotlib.gridspec.GridSpecBase` rows and columns, setting
    `wspace` and `hspace` to ``0``, and masking out every other row/column
    of the `~matplotlib.gridspec.GridSpecBase`, so they act as "spaces".
    These "spaces" are allowed to vary in width using the builtin
    `width_ratios` and `height_ratios` keyword args.
    """
    def __init__(self, nrows, ncols, **kwargs):
        """
        Parameters
        ----------
        nrows, ncols : int
            Number of rows, columns on the subplot grid.
        wspace, hspace : float or list of float
            The horizontal, vertical spacing between columns, rows of
            subplots. In `~proplot.subplots.subplots`, ``wspace``
            and ``hspace`` are in physical units. When calling
            `FlexibleGridSpecBase` directly, values are scaled relative to
            the average subplot width or height. If list, length of
            ``wspace``, ``hspace`` must be ``ncols-1``, ``nrows-1``.
        height_ratios, width_ratios : list of float
            Ratios for the width/height of columns/rows of subplots.
            For example, ``width_ratios=[1,2]`` specifes 2 columns of subplots,
            the second one twice as wide as the first.
        left, right, top, bottom : float or str
            Passed to `~matplotlib.gridspec.GridSpec`, denote the margin
            positions for the subplot grid in figure-relative coordinates.
        """
        self._nrows = nrows*2-1 # used with get_geometry
        self._ncols = ncols*2-1
        self._nrows_visible = nrows
        self._ncols_visible = ncols
        wratios, hratios, kwargs = self._spaces_as_ratios(**kwargs)
        super().__init__(self._nrows, self._ncols,
                hspace=0, wspace=0, # we implement these as invisible rows/columns
                width_ratios=wratios,
                height_ratios=hratios,
                **kwargs,
                )

    def __getitem__(self, key):
        """Magic obfuscation that renders `~matplotlib.gridspec.GridSpecBase`
        rows and columns designated as 'spaces' inaccessible."""
        # Get indices
        nrows, ncols = self.get_geometry()
        nrows_visible, ncols_visible = self.get_visible_geometry()
        if not isinstance(key, tuple): # usage gridspec[1,2]
            num1, num2 = _normalize(key, nrows_visible * ncols_visible)
        else:
            if len(key) == 2:
                k1, k2 = key
            else:
                raise ValueError(f'Invalid index {key!r}.')
            num1 = _normalize(k1, nrows_visible)
            num2 = _normalize(k2, ncols_visible)
            num1, num2 = np.ravel_multi_index((num1, num2), (nrows, ncols))
        # Correct for negative nums
        num1, num2 = _positem(num1), _positem(num2)
        return mgridspec.SubplotSpec(self, num1, num2)

    def _spaces_as_ratios(self,
        hspace=None, wspace=None, # spacing between axes
        height_ratios=None, width_ratios=None,
        **kwargs):
        """For keyword arg usage, see `FlexibleGridSpecBase`."""
        # Parse flexible input
        nrows, ncols = self.get_visible_geometry()
        hratios = np.atleast_1d(_notNone(height_ratios, 1))
        wratios = np.atleast_1d(_notNone(width_ratios,  1))
        hspace = np.atleast_1d(_notNone(hspace, np.mean(hratios)*0.10)) # this is relative to axes
        wspace = np.atleast_1d(_notNone(wspace, np.mean(wratios)*0.10))
        if len(wspace) == 1:
            wspace = np.repeat(wspace, (ncols-1,)) # note: may be length 0
        if len(hspace) == 1:
            hspace = np.repeat(hspace, (nrows-1,))
        if len(wratios) == 1:
            wratios = np.repeat(wratios, (ncols,))
        if len(hratios) == 1:
            hratios = np.repeat(hratios, (nrows,))

        # Verify input ratios and spacings
        # Translate height/width spacings, implement as extra columns/rows
        if len(hratios) != nrows:
            raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')
        if len(wratios) != ncols:
            raise ValueError(f'Got {ncols} columns, but {len(wratios)} wratios.')
        if len(wspace) != ncols-1:
            raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(wspace)}.')
        if len(hspace) != nrows-1:
            raise ValueError(f'Require {nrows-1} height spacings for {nrows} rows, got {len(hspace)}.')

        # Assign spacing as ratios
        # Also return extra kwargs, will be passed to superclass initializers
        # or superclass update method.
        nrows, ncols = self.get_geometry()
        wratios_final = [None]*ncols
        wratios_final[::2] = [*wratios]
        if ncols > 1:
            wratios_final[1::2] = [*wspace]
        hratios_final = [None]*nrows
        hratios_final[::2] = [*hratios]
        if nrows > 1:
            hratios_final[1::2] = [*hspace]
        return wratios_final, hratios_final, kwargs # bring extra kwargs back

    def get_wratios(self):
        """Alias for `~matplotlib.gridspec.Gridspecbase.get_width_ratios`."""
        return super().get_width_ratios()

    def get_hratios(self):
        """Alias for `~matplotlib.gridspec.Gridspecbase.get_height_ratios`."""
        return super().get_height_ratios()

    def get_visible_geometry(self):
        """Like `~matplotlib.gridspec.GridspecBase.get_geometry`, but returns
        the number of visible rows and columns, i.e. the number of rows and
        columns that aren't skipped over by
        `~FlexibleGridSpecBase.__getitem__`."""
        return self._nrows_visible, self._ncols_visible

class FlexibleGridSpec(FlexibleGridSpecBase, mgridspec.GridSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpec`."""
    def __init__(self, figure, **kwargs):
        """
        Parameters
        ----------
        figure : `Figure`
            The figure instance filled by this gridspec. Unlike
            `~matplotlib.gridspec.GridSpec`, this argument is required.
        **kwargs
            Passed to `~matplotlib.gridspec.GridSpec`.
        """
        super().__init__(figure=figure, **kwargs)

    def update(self, **kwargs):
        """Updates the width ratios, height ratios, and spacing for subplot
        columns and rows. The default `~matplotlib.gridspec.GridSpec.update`
        fails to update axes positions after successive edits to the figure,
        evidently because the figure is removed from the figure manager(s)."""
        # Convert spaces to ratios
        wratios, hratios, kwargs = self._spaces_as_ratios(**kwargs)
        self.set_width_ratios(wratios)
        self.set_height_ratios(hratios)
        # Validate args
        kwargs.pop('ncols', None)
        kwargs.pop('nrows', None)
        self.left   = kwargs.pop('left', None)
        self.right  = kwargs.pop('right', None)
        self.bottom = kwargs.pop('bottom', None)
        self.top    = kwargs.pop('top', None)
        if kwargs:
            raise ValueError(f'Unknown keyword arg(s): {kwargs}.')
        # Apply to figure and all axes
        fig = self.figure
        fig.subplotpars.update(self.left, self.bottom, self.right, self.top)
        for ax in fig.axes:
            ax.update_params()
            ax.set_position(ax.figbox)
        fig.stale = True

class FlexibleGridSpecFromSubplotSpec(FlexibleGridSpecBase, mgridspec.GridSpecFromSubplotSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpecFromSubplotSpec`."""
    pass

#-----------------------------------------------------------------------------#
# Helper funcs
#-----------------------------------------------------------------------------#
def _subplots_geometry(**kwargs):
    """Saves arguments passed to `subplots`, calculates gridspec settings and
    figure size necessary for requested geometry, and returns keyword args
    necessary to reconstruct and modify this configuration. Note that
    `wspace`, `hspace`, `left`, `right`, `top`, and `bottom` always have fixed
    physical units, then we scale figure width, figure height, and width
    and height ratios to accommodate requested geometry. And panel widths and
    heights are also fixed."""
    # Dimensions and geometry
    nrows, ncols       = kwargs['nrows'], kwargs['ncols']
    aspect, xref, yref = kwargs['aspect'], kwargs['xref'], kwargs['yref']
    width, height      = kwargs['width'], kwargs['height']
    axwidth, axheight  = kwargs['axwidth'], kwargs['axheight']
    # Gridspec settings
    wspace, hspace   = kwargs['wspace'], kwargs['hspace']
    wratios, hratios = kwargs['wratios'], kwargs['hratios']
    left, bottom     = kwargs['left'], kwargs['bottom']
    right, top       = kwargs['right'], kwargs['top']

    # Separate the panel and axes ratios
    axwspace, axhspace = wspace[3:-1:3], hspace[3:-1:3] # ignores figure panels
    axwidths, axheights = wratios[2:-1:3], hratios[2:-1:3]
    lpanels, tpanels = wratios[1:-1:3], hratios[1:-1:3]
    rpanels, bpanels = wratios[3:-1:3], hratios[3:-1:3]
    lpanel, rpanel = wratios[0], wratios[-1]
    tpanel, bpanel = hratios[0], hratios[-1]

    # Get reference properties, account for panel slots in space and ratios
    dx, dy = xref[1] - xref[0] + 1, yref[1] - yref[0] + 1
    rwspace = sum(axwspace[xref[0]:xref[1]])
    rhspace = sum(axhspace[yref[0]:yref[1]])
    rwratio = ncols*(sum(axwidths[xref[0]:xref[1]+1])/dx)/sum(axwidths)
    rhratio = nrows*(sum(axheights[yref[0]:yref[1]+1])/dy)/sum(axheights)
    if np.iterable(aspect):
        aspect = aspect[0]/aspect[1]

    # Determine figure and axes dims from input in width or height dimenion.
    # For e.g. common use case [[1,1,2,2],[0,3,3,0]], make sure we still scale
    # the reference axes like square even though takes two columns of gridspec!
    auto_width = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    if width is None and height is None: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = units(rc['subplots.axwidth'])
        if axheight is not None:
            auto_width = True
            axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows
            height = (axheight_all + top + bottom + sum(hspace)
                + sum(bpanels) + sum(tpanels) + bpanel + tpanel)
        if axwidth is not None:
            auto_height = True
            axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
            width = (axwidth_all + left + right + sum(wspace)
                + sum(lpanels) + sum(rpanels) + lpanel + rpanel)
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
    else:
        if height is not None:
            axheight_all = (height - top - bottom - sum(hspace)
                - sum(bpanels) - sum(tpanels) - bpanel - tpanel)
            axheight = (axheight_all*dy/nrows) + rhspace
        if width is not None:
            axwidth_all = (width - left - right - sum(wspace)
                - sum(lpanels) - sum(rpanels) - lpanel - rpanel)
            axwidth = (axwidth_all*dx/ncols) + rwspace

    # Automatically figure dim that was not specified above
    if auto_width:
        axwidth = axheight*aspect
        axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
        width = (axwidth_all + left + right + sum(wspace)
            + sum(lpanels) + sum(rpanels) + lpanel + rpanel)
    elif auto_height:
        axheight = axwidth/aspect
        axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows
        height = (axheight_all + top + bottom + sum(hspace)
            + sum(bpanels) + sum(tpanels) + bpanel + tpanel)
    if axwidth_all < 0:
        raise ValueError(f"Not enough room for axes (would have width {axwidth_all}). Try using tight=False, increasing figure width, or decreasing 'left', 'right', or 'wspace' spaces.")
    if axheight_all < 0:
        raise ValueError(f"Not enough room for axes (would have height {axheight_all}). Try using tight=False, increasing figure height, or decreasing 'top', 'bottom', or 'hspace' spaces.")

    # Reconstruct the ratios array with physical units for subplot slots
    # The panel slots are unchanged because panels have fixed widths
    axwidths = axwidth_all*np.array(axwidths)/sum(axwidths)
    axheights = axheight_all*np.array(axheights)/sum(axheights)
    wratios = np.repeat((None,), 3*ncols + 2)
    hratios = np.repeat((None,), 3*nrows + 2)
    wratios[2:-1:3], hratios[2:-1:3] = axwidths, axheights
    wratios[3:-1:3], hratios[3:-1:3] = rpanels, bpanels
    wratios[1:-1:3], hratios[1:-1:3] = lpanels, tpanels
    wratios[0], wratios[-1] = lpanel, rpanel
    hratios[0], hratios[-1] = tpanel, bpanel

    # Convert margins to figure-relative coordinates
    left, right = left/width,     1 - right/width
    top, bottom = 1 - top/height, bottom/height

    # Return gridspec keyword args
    gridspec_kw = {
        'ncols': ncols*3 + 2, 'nrows': nrows*3 + 2,
        'wspace': wspace, 'hspace': hspace,
        'width_ratios': wratios, 'height_ratios': hratios,
        'left': left, 'bottom': bottom, 'right': right, 'top': top,
        }

    return (width, height), gridspec_kw, kwargs

def _panels_kwargs(sides, kwargs, figure=False):
    """Converts global keywords like `space` and `width` to side-local
    keywords like `lspace` and `lwidth`, and applies default settings."""
    # Detect *unknown* keyword args, raise error if found!
    if not {*sides} <= {*'lrbt'}:
        raise ValueError(f'Invalid panel spec {sides!r}. Valid characters are "l", "r", "b", "t".')
    unknown = {}
    names = ('stack', 'share', 'space', 'width', 'sep')
    names = (('array', *names) if figure else names)
    regex = re.compile(f'^[tlrb]?({"|".join(names)})$')
    for key in (*kwargs.keys(),):
        if not regex.match(key):
            unknown[key] = kwargs.pop(key)
    if unknown:
        raise ValueError(f'Unknown keyword arg(s): {unknown}.')

    # Detect *ignored* panel keyword args, issue warning
    offsides = ''.join({*'lrbt'} - {*sides})
    regex = re.compile(f'^[{offsides}]({"|".join(names)})$')
    ikw = {key:value for key,value in kwargs.items() if regex.match(offsides)}
    if ikw:
        warnings.warn(f'Ignoring keyword arg(s): {ikw}. Active sides are: {sides!r}.')

    # Get default keyword values, add user input to orig dictionary
    kwout, kworig = {}, {}
    for side in sides:
        # Get arguments, try to use local first, global if no local found
        get = lambda key: _notNone(kwargs.get(side + key, None), kwargs.get(key, None))
        width = np.atleast_1d(units(get('width')))
        share = _notNone(get('share'), True)
        stack = _notNone(get('stack'), len(width))
        space = units(get('space'))
        sep = np.atleast_1d(units(get('sep')))

        # Validate stackable args, then store originals!
        if stack < 1:
            raise ValueError(f'{side+stack!r} argument must be integer >=1.')
        if len(width) == 1:
            width = np.repeat(width, (stack,))
        if len(width) != stack:
            raise ValueError(f'For side {side!r}, have {stack} stacked panels, but got {len(width)} widths.')
        if len(sep) == 1:
            sep = np.repeat(sep, (stack-1,))
        if len(sep) != stack-1:
            raise ValueError(f'For side {side!r}, have {stack} stacked panels, but got {len(sep)} separations.')
        kworig[side + 'space'] = space
        kworig[side + 'width'] = width
        kworig[side + 'sep'] = sep

        # Get default arg values
        space = _notNone(space, units(rc['subplots.' + ('panel' if share
            and not figure
            else 'xlab' if side == 'b' else 'ylab' if side == 'l'
            else 'inner' if figure else 'panel') + 'space']))
        width = np.array(width) # make a copy
        width[width==None] = units(rc['subplots.panelwidth'])
        sep = np.array(sep) # make a copy
        sep[sep==None] = units(rc['subplots.' + ('inner' if share
            else 'ylab' if side in 'lr' else 'xlab') + 'space'])

        # Store in output dictionary
        kwout[side + 'share'] = share
        kwout[side + 'width'] = width
        kwout[side + 'space'] = space
        kwout[side + 'sep'] = sep
        if figure:
            kwout[side + 'array'] = get('array')

    return kwout, kworig

#-----------------------------------------------------------------------------#
# Figure class and helper classes
#-----------------------------------------------------------------------------#
class _unlocker(object):
    """Suppresses warning message when adding subplots, and cleanly resets
    lock setting if exception raised."""
    def __init__(self, fig):
        self._fig = fig
    def __enter__(self):
        self._fig._locked = False
    def __exit__(self, *args):
        self._fig._locked = True

class _hidelabels(object):
    """Hides objects temporarily so they are ignored by the tight bounding box
    algorithm."""
    def __init__(self, *args):
        self._labels = args
    def __enter__(self):
        for label in self._labels:
            label.set_visible(False)
    def __exit__(self, *args):
        for label in self._labels:
            label.set_visible(True)

class Figure(mfigure.Figure):
    """The `~matplotlib.figure.Figure` class returned by `subplots`. At
    draw-time, an improved "tight layout" adjustment is applied, and
    the space around the figure edge, between subplots, and between
    panels is changed to accommodate subplot content. Figure dimensions
    may be automatically scaled to preserve subplot aspect ratios."""
    def __init__(self,
        tight=None,
        pad=None, axpad=None, panelpad=None,
        autoformat=True, ref=1,
        subplots_kw=None, gridspec_kw=None, subplots_orig_kw=None,
        tight_layout=None, constrained_layout=None,
        **kwargs):
        """
        Parameters
        ----------
        tight : bool, optional
            Toggles automatic tight layout adjustments.
        pad : float or str, optional
            Padding around edge of figure. Defaults to ``rc['subplots.pad']``.
            If float, units are inches. If string, units are interpreted by
            `~proplot.utils.units`.
        axpad : float or str, optional
            Padding between subplots in adjacent columns and rows. Defaults to
            ``rc['subplots.axpad']``. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        panelpad : float or str, optional
            Padding between subplots and axes panels, and between "stacked"
            panels. If float, units are inches. If string, units are
            interpreted by `~proplot.utils.units`.
        autoformat : bool, optional
            Whether to automatically format the axes when a `~pandas.Series`,
            `~pandas.DataFrame` or `~xarray.DataArray` is passed to a plotting
            command.
        gridspec_kw, subplots_kw, subplots_orig_kw
            Keywords used for initializing the main gridspec, for initializing
            the figure, and original spacing keyword args used for initializing
            the figure that override tight layout spacing.

        Other parameters
        ----------------
        tight_layout, constrained_layout
            Ignored, because ProPlot uses its own tight layout algorithm.
        **kwargs
            Passed to `matplotlib.figure.Figure`.
        """
        # Initialize first, because need to provide fully initialized figure
        # as argument to gridspec, because matplotlib tight_layout does that
        if tight_layout or constrained_layout:
            warnings.warn(f'Ignoring tight_layout={tight_layout} and contrained_layout={constrained_layout}. ProPlot uses its own tight layout algorithm, activated by default or with tight=True.')
        super().__init__(**kwargs)
        self._locked = False
        self._pad = units(_notNone(pad, rc['subplots.pad']))
        self._axpad = units(_notNone(axpad, rc['subplots.axpad']))
        self._panelpad = units(_notNone(panelpad, rc['subplots.panelpad']))
        self._auto_format = autoformat
        self._auto_tight_layout = _notNone(tight, rc['tight'])
        self._ref_num = ref
        self._main_axes = []
        self._subplots_orig_kw = subplots_orig_kw
        self._subplots_kw = subplots_kw
        self._leftpanel   = axes.EmptyPanel()
        self._bottompanel = axes.EmptyPanel()
        self._rightpanel  = axes.EmptyPanel()
        self._toppanel    = axes.EmptyPanel()
        self._main_gridspec = FlexibleGridSpec(self, **(gridspec_kw or {'nrows':1, 'ncols':1}))
        self.suptitle('') # add _suptitle attribute

    def _align_adjust(self, renderer):
        """Aligns row labels, column labels, and super titles, and applies
        tight layout automatic spacing."""
        self._adjust_aspect()
        self._align_axislabels(False)
        self._align_suplabels(renderer)
        if self._auto_tight_layout:
            self._adjust_tight_layout(renderer)
        self._align_axislabels(True)

    def _adjust_aspect(self):
        """Adjust average aspect ratio used for gridspec calculations. This
        fixes grids with identically fixed aspect ratios, e.g. identically
        zoomed-in cartopy projections and imshow images."""
        # Get aspect ratio
        if not self._main_axes:
            return
        ax = self._main_axes[self._ref_num-1]
        mode = ax.get_aspect()
        aspect = None
        if mode == 'equal':
            xscale, yscale = ax.get_xscale(), ax.get_yscale()
            if xscale == 'linear' and yscale == 'linear':
                aspect = 1.0/ax.get_data_ratio()
            elif xscale == 'log' and yscale == 'log':
                aspect = 1.0/ax.get_data_ratio_log()
            else:
                pass # matplotlib issues warning, forces aspect == 'auto'
        # Apply aspect
        # Account for floating point errors
        if aspect is not None:
            aspect = round(aspect*1e10)*1e-10
            subplots_kw = self._subplots_kw
            aspect_prev = round(subplots_kw['aspect']*1e10)*1e-10
            if aspect != aspect_prev:
                subplots_kw['aspect'] = aspect
                self._subplots_geometry()

    @_counter
    def _adjust_tight_layout(self, renderer=None,
        borders=True, subplots=True, panels=True
        ):
        """Applies tight layout scaling that permits flexible figure
        dimensions and preserves panel widths and subplot aspect ratios.
        The `renderer` should be a `~matplotlib.backend_bases.RendererBase`
        instance. If ``None``, renderer is inferred from
        `~matplotlib.figure.Figure.canvas.get_renderer`."""
        # Initial stuff
        axs = self._iter_axes()
        gridspec = self._main_gridspec
        subplots_kw = self._subplots_kw
        subplots_orig_kw = self._subplots_orig_kw # tight layout overrides
        if not axs or not subplots_kw or not subplots_orig_kw:
            return

        # Tight box *around* figure
        # Get bounds from old bounding box
        pad = self._pad
        obox = self.bbox_inches # original bbox
        bbox = self.get_tightbbox(renderer)
        l, b = bbox.xmin, bbox.ymin # left bottom margins
        r, t = obox.xmax - bbox.xmax, obox.ymax - bbox.ymax # top right margin *deltas*
        # Apply new bounds, permitting user overrides
        # TODO: Account for bounding box NaNs?
        for key,offset in zip(('left','right','top','bottom'),(l,r,t,b)):
            previous = subplots_orig_kw[key]
            current = subplots_kw[key]
            subplots_kw[key] = _notNone(previous, current - offset + pad)

        # Get arrays storing gridspec spacing args
        axpad = self._axpad
        panelpad = self._panelpad
        nrows, ncols = gridspec.get_visible_geometry()
        wspace, hspace = subplots_kw['wspace'], subplots_kw['hspace']
        wspace_orig = subplots_orig_kw['wspace']
        hspace_orig = subplots_orig_kw['hspace']
        wseps_orig = subplots_orig_kw['wseps']
        hseps_orig = subplots_orig_kw['hseps']
        # Get new subplot spacings, axes panel spacing, figure panel spacing
        spaces, ratios = [], [] # update these
        wratios = np.array(gridspec.get_wratios()[::2])
        hratios = np.array(gridspec.get_hratios()[::2])
        for (wh, xy, ixy, nacross,
            iratios, ispace,
            ispace_orig, iseps_orig) in zip('wh', 'xy', 'yx',
            (nrows,ncols), (wratios,hratios), (wspace,hspace),
            (wspace_orig,hspace_orig), (wseps_orig,hseps_orig),
            ):
            # Iterate along spaces
            jspace = [*ispace]
            ralong = np.array([ax._range_gridspec(xy, True) for ax in axs])
            racross = np.array([ax._range_gridspec(ixy, True) for ax in axs])
            for i,(space,space_orig) in enumerate(zip(ispace,ispace_orig)):
                # Special case where (1) current column is zero-width, so want
                # zero space next to it, or (2) no more non-zero columns
                # to the right of the current column
                pad = axpad
                if i == 0: # left/top figure panel
                    off1 = 0
                    off2 = (1 if iratios[i+1] == 0 else 0)
                    if iratios[i] == 0:
                        continue
                elif i == (len(ispace) - 1): # bottom/right figure panel
                    off2 = 0
                    off1 = (-1 if iratios[i] == 0 else 0)
                    if iratios[i+1] == 0:
                        continue
                elif (i % 3 != 0):
                    pad = panelpad # axes panels
                    off1, off2 = 0, 0
                    if iratios[i] == 0 or iratios[i+1] == 0: # column on either side of interface
                        continue
                else: # may have to jump across multiple columns
                    off1 = (-1 if iratios[i] == 0 else 0)
                    off2 = (1 if iratios[i+1] == 0 else 0)
                # Find axes that abutt aginst this space on each row
                # TODO: Fix offsets!!!
                groups = []
                filt1 = ralong[:,1] == i + off1 # i.e. right/bottom edge abutts against this space
                filt2 = ralong[:,0] == i + off2 + 1 # i.e. left/top edge abutts against this space
                for j in range(nacross): # e.g. each row
                    filt = ((racross[:,0] <= j) & (j <= racross[:,1]))
                    if sum(filt) < 2: # no interface here
                        continue
                    idx1, = np.where(filt & filt1)
                    idx2, = np.where(filt & filt2)
                    if not idx1.size or not idx2.size:
                        continue
                    # Put these axes into unique groups. Store groups as
                    # (left axes, right axes) or (bottom axes, top axes) pairs.
                    # NOTE: If have stacked panels, there may be more than one
                    # axes in the same gridspec slot! Account for multiple matches
                    axs1, axs2 = [axs[k] for k in idx1], [axs[k] for k in idx2]
                    if xy != 'x':
                        axs1, axs2 = axs2, axs1 # yrange is top-to-bottom, so make this bottom-to-top
                    newgroup = True
                    for (group1,group2) in groups:
                        if (any(ax1 in group1 for ax1 in axs1) or
                            any(ax2 in group2 for ax2 in axs2)):
                            newgroup = False
                            group1.update(axs1)
                            group2.update(axs2)
                            break
                    if newgroup:
                        groups.append([{*axs1}, {*axs2}]) # form new group
                # Get spaces
                # Remember layout is lspace, lspaces[0], rspaces[0], wspace, ...
                # so panels spaces are located where i % 3 is 1 or 2
                jspaces = []
                for (group1,group2) in groups:
                    x1 = max(ax._range_tightbbox(xy)[1] for ax in group1)
                    x2 = min(ax._range_tightbbox(xy)[0] for ax in group2)
                    jspaces.append((x2 - x1)/self.dpi)
                if jspaces:
                    space = max(0, space - min(jspaces) + pad) # TODO: why max 0?
                    space = _notNone(space_orig, space) # only if user did not provide original space!!!
                jspace[i] = space
            spaces.append(jspace)

            # Next handle separation in stacked panel gridspecs
            # We are iterating through *ratios* this time
            # TODO: This is clunkier than before, maybe go back to iterating
            # over axes instead of over the space/ratios arrays?
            jratios = [*iratios]
            for i,(ratio,sep_orig) in enumerate(zip(iratios,iseps_orig)):
                # If this is not a panel column, no panels are activated, sep
                # is [] (i.e. no stack), or sep provided manually, continue
                panel = (i % 3 != 2) or i in (0,len(iratios)-1)
                if not panel or ratio == 0 or all(sep_orig!=None):
                    continue
                # Next
                idx, = np.where(ralong[:,0] == i)
                iaxs = [axs[k] for k in idx]
                if not iaxs:
                    warnings.warn(f'No panels found in column or row {i}. This is weird.')
                    continue
                # Check for geometry inconsistencies
                gridspecs = [iax.get_subplotspec().get_gridspec() for iax in iaxs]
                geometry = {tuple(gridspec.get_visible_geometry()) for gridspec in gridspecs}
                if len(geometry)>1:
                    warnings.warn(f'Differing numbers of stacked panels in column or row, skipping tight layout adjustments.')
                    continue
                # Get indices relative to GridSpecFromSubplotSpec
                # Also make sure we ahve *active* axes in each column of stack
                # for each row, each row of stack for each column
                idxs = [iax._range_gridspec(xy, False)[0] for iax in iaxs]
                counts = [idxs.count(idx) for idx in {*idxs}]
                if any(count!=counts[0] for count in counts):
                    warnings.warn(f'Differing numbers of stacked panels in column or row, skipping tight layout adjustments.')
                    continue
                # Finally get tight layout, figure out separation between
                # each column in stack
                pad = panelpad
                iadd = []
                idxs_sorted = sorted({*idxs})[:-1]
                for idx_i in idxs_sorted:
                    # Calculate necessary change to seps; ccounts for fact that
                    # gridspec coords are top-to-bottom, bounding box opposite
                    if sep_orig[idx_i] is not None: # TODO: check!
                        iadd.append(0)
                    else:
                        off1, off2 = (0, 1) if xy == 'x' else (1, 0)
                        x1 = [iax._range_tightbbox(xy)[1] for idx,iax in
                                zip(idxs,iaxs) if idx == idx_i + off1]
                        x2 = [iax._range_tightbbox(xy)[0] for idx,iax in
                                zip(idxs,iaxs) if idx == idx_i + off2]
                        iadd.append(pad - (min(x2) - max(x1))/self.dpi)
                # Adjust spacing by adjusting subplot and inner ratios arrays
                # NOTE: For inner ratios adjustment, we apply change to mutable
                # ratios list, and change is reflected in next gridspec update
                jratios[i] += sum(iadd)
                for idx,iax in zip(idxs,iaxs):
                    if idx not in idxs_sorted: # i.e. is on bottom/left
                        continue
                    igridspec = iax.get_subplotspec().get_gridspec()
                    kratios = getattr(igridspec, 'get_' + wh + 'ratios')()
                    kratios[idx*2+1] += iadd[idx]
            ratios.append(jratios)

        # Apply new spaces and widths
        subplots_kw.update({
            'wspace':spaces[0], 'hspace':spaces[1],
            'wratios':ratios[0], 'hratios':ratios[1]
            })
        self._subplots_geometry()

    @_counter
    def _add_panel(self, ax, side, order='C', **kwargs):
        """Hidden method that powers `~proplot.axes.panel_axes`. Makes more sense
        to define in subplots, because it does a bunch of alignment steps
        that rely on the figure instance."""
        # Checks
        s = side[0]
        side = _side_translate.get(s, s)
        if s not in 'lrbt':
            raise ValueError(f'Invalid side {side!r}.')
        gridspec = self._main_gridspec
        if gridspec is None:
            raise ValueError(f'Gridspec attribute is empty.')
        paxs = getattr(ax, s + 'panel')
        if paxs:
            warnings.warn(f'{side}panel already exists.')
            return paxs # return existing panels! should already be synced
        # Parse keyword args
        subplots_kw = self._subplots_kw
        subplots_orig_kw = self._subplots_orig_kw
        kwout, kworig = _panels_kwargs(s, kwargs, figure=False)
        width_orig = kworig[s + 'width']
        space_orig = kworig[s + 'space']
        sep_orig = kworig[s + 'sep']
        width = kwout[s + 'width']
        share = kwout[s + 'share']
        space = kwout[s + 'space']
        sep = kwout[s + 'sep']
        # Check that this matches geometry of other panels in row or col
        axs = ax._get_side_axes(s)
        for iax in axs:
            ipaxs = getattr(iax, s + 'panel')
            if not ipaxs:
                continue
            nstack, instack = len(width), len(ipaxs)
            if nstack != instack:
                raise ValueError(f'You are trying to add {nstack} stacked panels to axes {iax}, but another axes in the same row/column has {instack} stacked panels.')

        # Find location in ratios and space arrays that we need to
        # modify to make room for the panel
        wh = ('w' if s in 'lr' else 'h')
        (x0, x1) = ax._range_gridspec('x', True)
        (y0, y1) = ax._range_gridspec('y', True)
        if s == 'l':
            x, y = x0-1, slice(y0,y1+1)
        elif s == 't':
            x, y = slice(x0,x1+1), y0-1
        elif s == 'r':
            x, y = x1+1, slice(y0,y1+1)
        else:
            x, y = slice(x0,x1+1), y1+1
        iratio = (x if s in 'lr' else y) # index in ratio array
        ispace = (iratio if s in 'lt' else iratio - 1) # index in space array

        # Default to original input space
        spaces_orig = subplots_orig_kw[wh + 'space']
        if spaces_orig[ispace] is None:
            spaces_orig[ispace] = space_orig
        else: # overwrite
            space = spaces_orig[ispace]
        # Default to original input width(s)
        widths_orig = subplots_orig_kw[wh + 'widths']
        iwidths_orig = widths_orig[iratio]
        if iwidths_orig is None:
            widths_orig[iratio] = width_orig
        else:
            mask = (width_orig != None)
            iwidths_orig[mask] = width_orig[mask]
            mask = (iwidths_orig != None)
            width[mask] = iwidths_orig[mask] # overwrite defaults with orig
        # Default to original input sep(s)
        seps_orig = subplots_orig_kw[wh + 'seps']
        iseps_orig = seps_orig[iratio]
        if iseps_orig is None:
            seps_orig[iratio] = sep_orig
        else:
            mask = (sep_orig != None)
            iseps_orig[mask] = sep_orig[mask]
            mask = (iseps_orig != None)
            sep[mask] = iseps_orig[mask] # overwrite

        # Apply to actual dict and get geometry
        subplots_kw[wh + 'ratios'][iratio] = sum(width) + sum(sep)
        subplots_kw[wh + 'space'][ispace] = space
        self._subplots_geometry()

        # Get subplotspec from main gridspec
        # NOTE: Critical that ratios are in physical units; we modify the
        # widths later on
        if s in 'lr':
            wspace, hspace = sep, []
            wratios, hratios = width, [1]
            nrows, ncols = 1, len(width)
        else:
            wspace, hspace = [], sep
            wratios, hratios = [1], width
            nrows, ncols = len(width), 1
        subplotspec = gridspec[y,x]
        # Put stack gridspec in subplotspec
        paxs = []
        igridspec = FlexibleGridSpecFromSubplotSpec(
            subplot_spec=subplotspec,
            nrows=nrows, ncols=ncols,
            wspace=wspace, hspace=hspace,
            width_ratios=wratios, height_ratios=hratios)
        # Draw panel(s)
        for i in range(len(width)):
            with self._unlock():
                pax = self.add_subplot(igridspec[i],
                    projection='panel',
                    side=side, parent=ax, share=share,
                    )
            paxs += [pax]

        # Add as axes_grid. Support 2D indexing, even though these are
        # always vector stacks, because consistency. See axes_grid docs.
        n = (1 if (s in 'tb' and order == 'C')
               or (s in 'lr' and order != 'C') else len(width))
        paxs = axes_grid(paxs, n=n, order=order)
        setattr(ax, '_' + side + 'panel', paxs) # hidden attribute for property
        # Share panel with axes from the main subplot, and update shared
        # axes settings externally
        if share:
            ax._share_panels_setup()
        self._share_axes_setup(ax)
        return paxs

    def _align_helper(self, xy, axs):
        """Gets figure coordinates for spanning labels and super title. The
        `xy` can be ``'x'`` or ``'y'``."""
        # Get position in figure relative coordinates
        ranges = np.array([ax._range_gridspec(xy, True) for ax in axs])
        min_, max_ = ranges[:,0].min(), ranges[:,1].max()
        axlo = axs[np.where(ranges[:,0] == min_)[0][0]]
        axhi = axs[np.where(ranges[:,1] == max_)[0][0]]
        lobox = axlo.get_subplotspec().get_position(self)
        hibox = axhi.get_subplotspec().get_position(self)
        if xy == 'x':
            pos = (lobox.x0 + hibox.x1)/2
        else:
            pos = (lobox.y1 + hibox.y0)/2 # 'lo' is actually on top, highest up in gridspec
        # Return axis suitable for spanning position
        spanax = axs[(np.argmin(ranges[:,0]) + np.argmax(ranges[:,1]))//2]
        return pos, spanax

    def _align_axislabels(self, b=True):
        """Aligns spanning *x* and *y* axis labels, accounting for figure
        margins and axes and figure panels."""
        # TODO: Ensure this is robust to complex panels and shared axes
        # NOTE: Need to turn off aligned labels before _adjust_tight_layout
        # call, so cannot put this inside Axes draw
        tracker = {*()}
        for ax in self._main_axes:
            if not isinstance(ax, axes.CartesianAxes):
                continue
            for xy,axis in zip('xy', (ax.xaxis, ax.yaxis)):
                s = axis.get_label_position()[0] # top or bottom, left or right
                span = getattr(ax, '_span' + xy)
                align = getattr(ax, '_align' + xy)
                if s not in 'bl' or axis in tracker:
                    continue
                axs = ax._get_side_axes(s)
                for _ in range(2):
                    axs = [getattr(ax, '_share' + xy) or ax for ax in axs]
                # Align axis label offsets
                axiss = [getattr(ax, xy + 'axis') for ax in axs]
                tracker.update(axiss)
                if span or align:
                    grp = getattr(self, '_align_' + xy + 'label_grp', None)
                    if grp is not None:
                        for ax in axs[1:]:
                            grp.join(axs[0], ax) # copied from source code, add to grouper
                    elif align:
                        warnings.warn(f'Aligning *x* and *y* axis labels required matplotlib >=3.1.0')
                if not span:
                    continue
                # Get spanning label position
                pos, spanax = self._align_helper(xy, axs)
                spanaxis = getattr(spanax, xy + 'axis')
                spanlabel = spanaxis.label
                if not hasattr(spanlabel, '_orig_transform'):
                    spanlabel._orig_transform = spanlabel.get_transform()
                    spanlabel._orig_position = spanlabel.get_position()
                if not b: # toggle off, done before tight layout
                    spanlabel.set_transform(spanlabel._orig_transform)
                    spanlabel.set_position(spanlabel._orig_position)
                    for axis in axiss:
                        axis.label.set_visible(True)
                else: # toggle on, done after tight layout
                    if xy == 'x':
                        position = (pos, 1)
                        transform = mtransforms.blended_transform_factory(
                                self.transFigure, mtransforms.IdentityTransform())
                    else:
                        position = (1, pos)
                        transform = mtransforms.blended_transform_factory(
                                mtransforms.IdentityTransform(), self.transFigure)
                    for axis in axiss:
                        axis.label.set_visible((axis is spanaxis))
                    spanlabel.update({'position':position, 'transform':transform})

    def _align_suplabels(self, renderer):
        """Adjusts position of row and column labels, and aligns figure
        super title accounting for figure marins and axes and figure panels."""
        # Offset using tight bounding boxes
        # TODO: Super labels fail with popup backend!! Fix this
        # NOTE: Must use get_tightbbox so (1) this will work if tight layout
        # mode if off and (2) actually need *two* tight bounding boxes when
        # labels are present: 1 not including the labels, used to position
        # them, and 1 including the labels, used to determine figure borders
        suptitle = self._suptitle
        suptitle_on = suptitle.get_text().strip()
        width, height = self.get_size_inches()
        for s in 'lrbt':
            # Geometry
            nrows, ncols = self._main_gridspec.get_visible_geometry()
            xy = ('x' if s in 'lr' else 'y')
            idx = (0 if s in 'lt' else 1)
            edge = (2 if s in 'lt' else ncols-3 if s == 'r' else nrows-3)
            # Get axes and offset the label to relevant panel
            axs = [ax._reassign_suplabel(s) for ax in self._main_axes
                   if ax._range_gridspec(xy, True)[idx] == edge]
            labels = [getattr(ax, '_' + s + 'label') for ax in axs]
            coords = [None]*len(axs)
            if s == 't' and suptitle_on:
                supaxs = axs
            with _hidelabels(*labels):
                for i,(ax,label) in enumerate(zip(axs,labels)):
                    label_on = label.get_text().strip()
                    if not label_on:
                        continue
                    # Get coord from tight bounding box
                    # Include twin axes and panels along the same side
                    extra = ('bt' if s in 'lr' else 'lr')
                    icoords = []
                    for iax in ax._iter_twins_panels(extra):
                        bbox = iax.get_tightbbox(renderer)
                        if s == 'l':
                            jcoords = (bbox.xmin, 0)
                        elif s == 'r':
                            jcoords = (bbox.xmax, 0)
                        elif s == 't':
                            jcoords = (0, bbox.ymax)
                        else:
                            jcoords = (0, bbox.ymin)
                        c = self.transFigure.inverted().transform(jcoords)
                        c = (c[0] if s in 'lr' else c[1])
                        icoords.append(c)
                    # Offset, and offset a bit extra for left/right labels
                    # See: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
                    fontsize = label.get_fontsize()
                    if s in 'lr':
                        scale1, scale2 = 0.6, width
                    else:
                        scale1, scale2 = 0.3, height
                    if s in 'lb':
                        coords[i] = min(icoords) - (scale1*fontsize/72)/scale2
                    else:
                        coords[i] = max(icoords) + (scale1*fontsize/72)/scale2
                # Assign coords
                coords = [i for i in coords if i is not None]
                if coords:
                    key = ('x' if s in 'lr' else 'y')
                    if s in 'lb':
                        c = min(coords)
                    else:
                        c = max(coords)
                    for label in labels:
                        label.update({key: c})

        # Update super title position
        if suptitle_on:
            ys = []
            for ax in supaxs:
                bbox = ax.get_tightbbox(renderer)
                _, y = self.transFigure.inverted().transform((0, bbox.ymax))
                ys.append(y)
            x, _ = self._align_helper('x', axs)
            y = max(ys) + (0.3*suptitle.get_fontsize()/72)/height
            kw = {'x':x, 'y':y, 'ha':'center', 'va':'bottom',
                  'transform':self.transFigure}
            suptitle.update(kw)

    def _unlock(self):
        """Prevents warning message when adding subplots one-by-one, used
        internally."""
        return _unlocker(self)

    def _update_axislabels(self, axis=None, **kwargs):
        """Applies axis labels to the relevant shared axis. If spanning
        labels are toggled, keeps the labels synced for all subplots in the
        same row or column. Label positions will be adjusted at draw-time
        with _align_axislabels."""
        # Sync text for axes encompassed by spanning labels
        # TODO: Make axis label alignment more robust when panels present
        xy = axis.axis_name
        if xy not in 'xy':
            return
        ax = axis.axes
        ax = getattr(ax, '_share' + xy, None) or ax # try to use panel axes
        s = axis.get_label_position()[0]  # top or bottom, left or right
        if s in 'lb' and getattr(ax, '_span' + xy):
            axs = ax._get_side_axes(s)
        else:
            axs = [ax]
        for ax in axs:
            for _ in range(2): # try up to two levels deep
                ax = getattr(ax, '_share' + xy, None) or ax
            axis = getattr(ax, xy + 'axis')
            axis.label.set_visible(True)
            axis.label.update(kwargs)

    def _update_suplabels(self, ax, side, labels, **kwargs):
        """Assigns side labels, updates label settings."""
        s = side[0]
        if s not in 'lrbt':
            raise ValueError(f'Invalid label side {side!r}.')
        if not self._main_gridspec:
            raise ValueError(f'Gridspec missing, cannot update labels.')
        # Geometry
        nrows, ncols = self._main_gridspec.get_visible_geometry()
        xy = ('x' if s in 'lr' else 'y')
        ixy = ('y' if s in 'lr' else 'x')
        idx = (0 if s in 'lt' else 1)
        edge = (2 if s in 'lt' else ncols-3 if s == 'r' else nrows-3)
        # Update label text for axes on the edge
        if not self._main_axes: # occurs if this is called while adding axes
            axs = [ax]
        else:
            axs = [ax for ax in self._main_axes if ax._range_gridspec(xy, True)[idx] == edge]
            axs = [ax for _,ax in sorted(zip([ax._range_gridspec(ixy, True)[0] for ax in axs], axs))] # order by yrange
        if labels is None or isinstance(labels, str): # common during testing
            labels = [labels]*len(axs)
        if len(labels) != len(axs):
            raise ValueError(f'Got {len(labels)} {s}labels, but there are {len(axs)} axes along that side.')
        for ax,label in zip(axs,labels):
            obj = getattr(ax, '_' + s + 'label')
            if label is not None and obj.get_text() != label:
                obj.set_text(label)
            if kwargs:
                obj.update(kwargs)

    def _update_suptitle(self, title, **kwargs):
        """Assign figure "super title"."""
        if title is not None and self._suptitle.get_text() != title:
            self._suptitle.set_text(title)
        if kwargs:
            self._suptitle.update(kwargs)

    def _share_axes_setup(self, ref=None):
        """Applies axis sharing to groups of axes that share the same
        horizontal or vertical extent."""
        axs = [ref] if ref is not None else self._main_axes
        # Share x axes
        ranges = {tuple(ax._range_gridspec('x', True)) for ax in axs}
        for irange in ranges:
            iaxs = [ax for ax in axs if tuple(ax._range_gridspec('x', True)) == irange]
            parent = iaxs.pop(np.argmax([iax._range_gridspec('y', True)[1] for iax in iaxs]))
            for child in iaxs:
                child._sharex_setup(parent, parent._sharex_level)
        # Share y axes
        ranges = {tuple(ax._range_gridspec('y', True)) for ax in axs}
        for irange in ranges:
            iaxs = [ax for ax in axs if tuple(ax._range_gridspec('y', True)) == irange]
            parent = iaxs.pop(np.argmin([iax._range_gridspec('x', True)[0] for iax in iaxs]))
            for child in iaxs:
                child._sharey_setup(parent, parent._sharey_level)

    def _subplots_geometry(self):
        """Sets up consistent geometry for subplots."""
        figsize, gridspec_kw, _ = _subplots_geometry(**self._subplots_kw)
        self.set_size_inches(figsize)
        self._main_gridspec.update(**gridspec_kw)

    def add_subplot(self, *args, **kwargs):
        """Issues warning for new users that try to call
        `~matplotlib.figure.Figure.add_subplot` or
        `~matplotlib.figure.Figure.colorbar` manually."""
        if self._locked:
            warnings.warn('Using "fig.add_subplot" or "fig.colorbar" with ProPlot figures may result in unexpected behavior. Please use "proplot.subplots" and "Axes.colorbar" instead.')
        ax = super().add_subplot(*args, **kwargs)
        return ax

    @_counter
    def draw(self, renderer):
        """Before drawing the figure, applies "tight layout" and aspect
        ratio-conserving adjustments, and aligns row and column labels."""
        # Renderer fixes
        # WARNING: *Critical* that draw() is invoked with the same renderer
        # FigureCanvasAgg.print_png() uses to render the image. But print_png()
        # calls get_renderer() after draw(), and get_renderer() returns a new
        # renderer if it detects renderer dims and figure dims are out of sync!
        # 1. Could use 'get_renderer' to update 'canvas.renderer' with the new
        #    figure width and height, then use that renderer for rest of draw
        #    This repair *breaks* just the *macosx* popup backend and not the
        #    qt backend! So for now just employ simple exception if this is
        #    macosx backend.
        # 2. Could set '_lastKey' on canvas and 'width' and 'height' on renderer,
        #    but then '_renderer' was initialized with wrong width and height,
        #    which causes bugs. And _renderer was generated with cython code
        # WARNING: Vector graphic renderers are another ballgame, *impossible*
        # to consistently apply successive figure size changes. SVGRenderer
        # and PDFRenderer both query the size in inches before calling draw,
        # and cannot modify PDFPage or SVG renderer props inplace, so idea was
        # to override get_size_inches. But when get_size_inches is called, the
        # canvas has no renderer, so cannot apply tight layout yet!
        self._align_adjust(renderer)
        canvas = getattr(self, 'canvas', None)
        if hasattr(canvas, 'get_renderer') and (mbackend is None or
            not isinstance(canvas, mbackend.FigureCanvasMac)):
            renderer = canvas.get_renderer()
            canvas.renderer = renderer
        return super().draw(renderer)

    def savefig(self, filename, **kwargs):
        """
        Before saving the figure, applies "tight layout" and aspect
        ratio-conserving adjustments, and aligns row and column labels.

        Parameters
        ----------
        filename : str
            The file path. User directories are automatically
            expanded, e.g. ``fig.save('~/plots/plot.png')``.
        **kwargs
            Passed to `~matplotlib.figure.Figure.savefig`.
        """
        filename = os.path.expanduser(filename)
        canvas = getattr(self, 'canvas', None)
        if hasattr(canvas, 'get_renderer'):
            renderer = canvas.get_renderer()
            canvas.renderer = renderer
            self._align_adjust(renderer)
        else:
            warnings.warn('Renderer unknown, could not adjust layout before saving.')
        super().savefig(filename, **kwargs)

    save = savefig
    """Alias for `~Figure.savefig`, because calling ``fig.savefig``
    is sort of redundant."""

    @property
    def leftpanel(self):
        """An `~proplot.subplots.axes_grid` of the left panel stack."""
        return self._leftpanel
    @property
    def rightpanel(self):
        """An `~proplot.subplots.axes_grid` of the right panel stack."""
        return self._rightpanel
    @property
    def toppanel(self):
        """An `~proplot.subplots.axes_grid` of the top panel stack."""
        return self._toppanel
    @property
    def bottompanel(self):
        """An `~proplot.subplots.axes_grid` of the bottom panel stack."""
        return self._bottompanel

    lpanel = leftpanel
    """Alias for `~Figure.leftpanel`."""
    rpanel = rightpanel
    """Alias for `~Figure.rightpanel`."""
    tpanel = toppanel
    """Alias for `~Figure.toppanel`."""
    bpanel = bottompanel
    """Alias for `~Figure.bottompanel`."""

    def _iter_axes(self):
        """Iterates over all axes and panels in the figure belonging to the
        `~proplot.axes.Axes` class. Excludes inset and twin axes."""
        axs = []
        for ax in (*self._main_axes, *self.leftpanel, *self.rightpanel,
                   *self.bottompanel, *self.toppanel):
            if not ax or not ax.get_visible():
                continue
            axs.append(ax)
        for ax in axs:
            for s in 'lrbt':
                for iax in getattr(ax, s + 'panel'):
                    if not iax or not iax.get_visible():
                        continue
                    axs.append(iax)
        return axs

#-----------------------------------------------------------------------------#
# Primary plotting function, used to create figure/axes
#-----------------------------------------------------------------------------#
def _journals(journal):
    """Journal sizes for figures."""
    table = {
        'pnas1': '8.7cm', # if 1 number specified, this is a tuple
        'pnas2': '11.4cm',
        'pnas3': '17.8cm',
        'ams1': 3.2, # spec is in inches
        'ams2': 4.5,
        'ams3': 5.5,
        'ams4': 6.5,
        'agu1': ('95mm',  '115mm'),
        'agu2': ('190mm', '115mm'),
        'agu3': ('95mm',  '230mm'),
        'agu4': ('190mm', '230mm'),
        }
    value = table.get(journal, None)
    if value is None:
        raise ValueError(f'Unknown journal figure size specifier {journal!r}. ' +
                          'Current options are: ' + ', '.join(table.keys()))
    # Return width, and optionally also the height
    width, height = None, None
    try:
        width, height = value
    except TypeError:
        width = value
    return width, height

def _axes_dict(naxs, value, kw=False, default=None):
    """Build a dictionary that looks like ``{1:value1, 2:value2, ...}`` or
    ``{1:{key1:value1, ...}, 2:{key2:value2, ...}, ...}`` for storing
    standardized axes-specific properties or keyword args."""
    # First build up dictionary
    # 1) 'string' or {1:'string1', (2,3):'string2'}
    if not kw:
        if np.iterable(value) and not isinstance(value, (str,dict)):
            value = {num+1: item for num,item in enumerate(value)}
        elif not isinstance(value, dict):
            value = {range(1, naxs+1): value}
    # 2) {'prop':value} or {1:{'prop':value1}, (2,3):{'prop':value2}}
    else:
        nested = [isinstance(value,dict) for value in value.values()]
        if not any(nested): # any([]) == False
            value = {range(1, naxs+1): value.copy()}
        elif not all(nested):
            raise ValueError('Pass either of dictionary of key value pairs or a dictionary of dictionaries of key value pairs.')
    # Then *unfurl* keys that contain multiple axes numbers, i.e. are meant
    # to indicate properties for multiple axes at once
    kwargs = {}
    for nums,item in value.items():
        nums = np.atleast_1d(nums)
        for num in nums.flat:
            if not kw:
                kwargs[num] = item
            else:
                kwargs[num] = item.copy()
    # Fill with default values
    for num in range(1, naxs+1):
        if num not in kwargs:
            if kw:
                kwargs[num] = {}
            else:
                kwargs[num] = default
    # Verify numbers
    if {*range(1, naxs+1)} != {*kwargs.keys()}:
        raise ValueError(f'Have {naxs} axes, but {value} has properties for axes {", ".join(str(i) for i in sorted(kwargs.keys()))}.')
    return kwargs

def subplots(array=None, ncols=1, nrows=1,
    ref=1, order='C',
    aspect=1, figsize=None,
    width=None,  height=None, axwidth=None, axheight=None, journal=None,
    wwidth=None, hwidth=None,
    axwidths=None, axheights=None,
    hspace=None, wspace=None, space=None,
    hratios=None, wratios=None,
    width_ratios=None, height_ratios=None,
    flush=None, wflush=None, hflush=None,
    left=None, bottom=None, right=None, top=None,
    tight=None, pad=None, axpad=None, panelpad=None,
    span=None, spanx=None, spany=None,
    align=None, alignx=None, aligny=None,
    share=None, sharex=None, sharey=None,
    panel=None, panels=None, axpanel=None, axpanels=None,
    axpanel_kw=None, axpanels_kw=None,
    basemap=False, proj=None, projection=None, proj_kw=None, projection_kw=None,
    autoformat=True,
    **kwargs):
    """
    Analogous to `matplotlib.pyplot.subplots`, creates a figure with a single
    axes or arbitrary grids of axes, any of which can be map projections,
    and optional "panels" along axes or figure edges.

    The parameters are sorted into the following rough sections: subplot grid
    specifications, figure and subplot sizes, axis sharing,
    figure panels, axes panels, and map projections.

    Parameters
    ----------
    ncols, nrows : int, optional
        Number of columns, rows. Ignored if `array` is not ``None``.
        Use these arguments for simpler subplot grids.
    order : {'C', 'F'}, optional
        Whether subplots are numbered in column-major (``'C'``) or row-major
        (``'F'``) order. Analogous to `numpy.array` ordering. This controls
        the order axes appear in the `axs` list, and the order of subplot
        a-b-c labeling (see `~proplot.axes.Axes.format`).
    array : array-like of int, optional
        2-dimensional array specifying complex grid of subplots. Think of
        this array as a "picture" of your figure. For example, the array
        ``[[1, 1], [2, 3]]`` creates one long subplot in the top row, two
        smaller subplots in the bottom row. Integers must range from 1 to the
        number of plots.

        ``0`` indicates an empty space. For example, ``[[1, 1, 1], [2, 0, 3]]``
        creates one long subplot in the top row with two subplots in the bottom
        row separated by a space.
    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    width, height : float or str, optional
        The figure width and height. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`. For example,
        ``width="10cm"`` creates a 10cm wide figure.
    journal : str, optional
        String name corresponding to an academic journal standard that is used
        to control the figure width (and height, if specified). Valid names
        are described in a table below.

    ref : int, optional
        The reference axes number. The `axwidth`, `axheight`, and `aspect`
        keyword args are applied to this axes, and aspect ratio is conserved
        for this axes in tight layout adjustment.
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
    hratios, wratios
        Aliases for `height_ratios`, `width_ratios`.
    width_ratios, height_ratios : float or list thereof, optional
        Passed to `FlexibleGridSpec`. The width
        and height ratios for the subplot grid. Length of `width_ratios`
        must match the number of rows, and length of `height_ratios` must
        match the number of columns.
    wspace, hspace, space : float or str or list thereof, optional
        Passed to `FlexibleGridSpec`, denotes the
        spacing between grid columns, rows, and both, respectively. If float
        or string, expanded into lists of length ``ncols-1`` (for `wspace`)
        or length ``nrows-1`` (for `hspace`). For each element of the list,
        if float, units are inches, and if string, units are interpreted by
        `~proplot.utils.units`.
    left, right, top, bottom : float or str, optional
        Passed to `FlexibleGridSpec`, denote the width of padding between the
        subplots and the figure edge. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

    sharex, sharey, share : {3, 2, 1, 0}, optional
        The "axis sharing level" for the *x* axis, *y* axis, or both
        axes. This can considerably redundancy in your figure.
        Options are as follows:

        0. No axis sharing. Also sets the default `spanx` and `spany` values
           to ``False``.
        1. Only draw *axis label* on the leftmost column (*y*) or
           bottommost row (*x*) of subplots. Axis tick labels
           still appear on every subplot.
        2. As in 1, but forces the axis limits to be identical. Axis
           tick labels still appear on every subplot.
        3. As in 2, but only show the *axis tick labels* on the
           leftmost column (*y*) or bottommost row (*x*) of subplots.

    spanx, spany, span : bool or {0, 1}, optional
        Default is ``False`` if `sharex`, `sharey`, or `share` are ``0``,
        ``True`` otherwise. Toggles "spanning" axis labels for the *x* axis,
        *y* axis, or both axes. When ``True``, a single, centered axis label
        is used for all axes with bottom and left edges in the same row or
        column.  This can considerably redundancy in your figure.

        "Spanning" labels integrate with "shared" axes. For example,
        for a 3-row, 3-column figure, with ``sharey > 1`` and ``spany=1``,
        your figure will have 1 ylabel instead of 9.
    alignx, aligny, align : bool or {0, 1}, optional
        Default is ``False``. Whether to `align axis labels
        <https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/align_labels_demo.html>`__
        for the *x* axis, *y* axis, or both axes. Only has an effect when
        `spanx`, `spany`, or `span` are ``False``.
    proj, projection : str or dict-like, optional
        The map projection name. The argument is interpreted as follows.

        * If string, this projection is used for all subplots. For valid
          names, see the :ref:`Table of projections`.
        * If list of string, these are the projections to use for each
          subplot in their `array` order.
        * If dict-like, keys are integers or tuple integers that indicate
          the projection to use for each subplot. If a key is not provided,
          that subplot will be a `~proplot.axes.CartesianAxes`. For example,
          in a 4-subplot figure, ``proj={2:'merc', (3,4):'stere'}``
          draws a Cartesian axes for the first subplot, a Mercator
          projection for the second subplot, and a Stereographic projection
          for the second and third subplots.

    proj_kw, projection_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` classes on instantiation.
        If dictionary of properties, applies globally. If *dictionary of
        dictionaries* of properties, applies to specific subplots, as
        with `proj`.

        For example, with ``ncols=2`` and
        ``proj_kw={1:{'lon_0':0}, 2:{'lon_0':180}}``, the projection in
        the left subplot is centered on the prime meridian, and the projection
        in the right subplot is centered on the international dateline.
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.
        If boolean, applies to all subplots. If dictionary, values apply to
        specific subplots, as with `proj`.

    panel : str, optional
        Adds :ref:`Global figure panels` to the specified side.
        String should contain any of the characters ``'l'`` (left), ``'r'``
        (right), ``'b'`` (bottom), or ``'t'`` (top). For example, ``'br'``
        will draw panels on the right and bottom sides of the figure.

        Figure panels are saved as the `~Figure.leftpanel`,
        `~Figure.rightpanel`, `~Figure.bottompanel`, or
        `~Figure.toppanel` attributes, with respective shorthands
        `~Figure.lpanel`, `~Figure.rpanel`, `~Figure.bpanel`, and
        `~Figure.tpanel`.
    panels : str, optional
        As with `panel`, but the default behavior is to assign a panel
        to *every* row or column of subplots. Individual panels can then
        be accessed with e.g. ``fig.leftpanel[0]``, ``fig.leftpanel[1]``.
    [lrbt]array : list of int, optional
        Defines how figure panels span rows and columns of subplots.
        Interpreted like `array` -- the integers specify panels that span
        *arbitrary, contiguous* columns or rows of subplots.

        For example, ``plot.suplots(ncols=3, panels='b', barray=[1, 2, 2])``
        draws a panel on the bottom of the first column and spanning the bottom
        of the right 2 columns, and ``barray=[0, 2, 2]`` only draws a panel
        underneath the right 2 columns -- as with `array`, the ``0`` indicates
        an empty space.
    [lrbt]width, [lrbt]share, [lrbt]stack, [lrbt]sep
        As in `~proplot.axes.Axes.panel_axes`. Use e.g. `lwidth`, `rwidth`,
        `bwidth`, and `twidth` to apply these settings to different panels.
    [lrtb]space : float or str, optional
        As in `~proplot.axes.Axes.panel_axes`, but controls space
        between the main subplot grid and the figure panels.

    Other parameters
    ----------------
    tight : bool, optional
        Toggles automatic tight layout adjustments for all spacings that you
        did not specify manually. For example, with ``left=0.1``, the right,
        top, bottom, and inner spaces will be determined automatically.
    pad, axpad, panelpad : float or str, optional
        Padding for automatic tight layout adjustments. See `Figure` for
        details.
    autoformat : bool, optional
        Whether to automatically format axes when special datasets are
        passed to plotting commands. See `Figure` for details.
    axpanel, axpanels : str or dict-like, optional
        Adds :ref:`Bulk axes panels` to subplots with
        `~proplot.axes.Axes.panel_axes`. Both `axpanel` and `axpanels`
        are acceptable. The argument is interpreted as follows.

        * If string, panels are drawn on the same side for all subplots.
          String should contain any of the characters ``'l'`` (left panel),
          ``'r'`` (right panel), ``'t'`` (top panel), or ``'b'`` (bottom
          panel). For example, ``'rt'`` will draw a right and top panel.
        * If dict-like, panels can be drawn on different sides for
          different subplots. For example, consider a 3-subplot figure.
          With ``axpanels={2:'r', 3:'l'}``, subplot 1 will have no panel,
          subplot 2 will have a panel on the right side, and subplot 3
          will have a panel on the left side.

    axpanel_kw, axpanels_kw : dict-like, optional
        Keyword args passed to `~proplot.axes.Axes.panel_axes` for panels
        listed in `axpanel` and `axpanels`. If dictionary of properties,
        applies globally. If *dictionary of dictionary* of properties, applies
        to specific subplots, as with `axpanels`.

        For example, consider a 2-subplot figure with ``axpanels='l'``.
        With ``axpanel_kw={'lwidth':1}``, both left panels will be 1 inch wide.
        With ``axpanel_kw={1:{'lwidth':1}, 2:{'lwidth':0.5}}``, the left
        subplot panel will be 1 inch wide and the right subplot panel will be
        0.5 inches wide.

    Returns
    -------
    f : `Figure`
        The figure instance.
    axs : `axes_grid`
        A special list of axes instances. See `axes_grid`.


    Current options for the `journal` keyword argument are as follows.
    If you'd like to add additional standards, feel free to submit a pull request

    ===========  ====================  ==========================================================================================================================================================
    Key          Size description      Organization
    ===========  ====================  ==========================================================================================================================================================
    ``'pnas1'``  1-column              `Proceedings of the National Academy of Sciences <http://www.pnas.org/page/authors/submission>`__
    ``'pnas2'``  2-column              ”
    ``'pnas3'``  landscape page        ”
    ``'ams1'``   1-column              `American Meteorological Society <https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/>`__
    ``'ams2'``   small 2-column        ”
    ``'ams3'``   medium 2-column       ”
    ``'ams4'``   full 2-column         ”
    ``'agu1'``   1-column              `American Geophysical Union <https://publications.agu.org/author-resource-center/figures-faq/>`__
    ``'agu2'``   2-column              ”
    ``'agu3'``   full height 1-column  ”
    ``'agu4'``   full height 2-column  ”
    ===========  ====================  ==========================================================================================================================================================
    """
    #-------------------------------------------------------------------------#
    # Array setup
    #-------------------------------------------------------------------------#
    rc._getitem_mode = 0 # ensure still zero; might be non-zero if had error in 'with context' block
    # Build array
    if order not in ('C','F'): # better error message
        raise ValueError(f'Invalid order {order!r}. Choose from "C" (row-major, default) and "F" (column-major).')
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        array = array.reshape((nrows, ncols), order=order)
    # Standardize array
    try:
        array = np.array(array, dtype=int) # enforce array type
        if array.ndim == 1:
            array = array[None,:] if order == 'C' else array[:,None] # interpret as single row or column
        elif array.ndim != 2:
            raise ValueError
        array[array == None] = 0 # use zero for placeholder
    except (TypeError,ValueError):
        raise ValueError(f'Invalid subplot array {array}. Must be 1D or 2D array of integers.')
    # Get other props
    nums = np.unique(array[array != 0])
    naxs = len(nums)
    if {*nums.flat} != {*range(1, naxs+1)}:
        raise ValueError('Invalid subplot array {array}. Numbers must span integers 1 to naxs (i.e. cannot skip over numbers), with 0 representing empty spaces.')
    if ref not in nums:
        raise ValueError(f'Invalid reference number {ref}. For array {array}, must be one of {nums}.')
    nrows, ncols = array.shape

    #-------------------------------------------------------------------------#
    # Shared and spanning axes, panel syncing
    #-------------------------------------------------------------------------#
    # Figure out rows and columns "spanned" by each axes in list, for
    # axis sharing and axis label spanning settings
    sharex = int(_notNone(sharex, share, rc['share']))
    sharey = int(_notNone(sharey, share, rc['share']))
    if sharex not in range(4) or sharey not in range(4):
        raise ValueError(f'Axis sharing level can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels). Got sharex={sharex} and sharey={sharey}.')
    spanx = _notNone(spanx, span, 0 if sharex == 0 else None, rc['span'])
    spany = _notNone(spany, span, 0 if sharey == 0 else None, rc['span'])
    alignx = _notNone(alignx, align)
    aligny = _notNone(aligny, align)
    if (spanx and alignx) or (spany and aligny):
        warnings.warn(f'The "alignx" and "aligny" args have no effect when "spanx" and "spany" are True.')
    alignx = _notNone(alignx, rc['align'])
    aligny = _notNone(alignx, rc['align'])
    # Get some axes properties, where locations are sorted by axes id.
    # NOTE: These ranges are endpoint exclusive, like a slice object!
    axids = [np.where(array == i) for i in np.sort(np.unique(array)) if i > 0] # 0 stands for empty
    xrange = np.array([[x.min(), x.max()] for _,x in axids])
    yrange = np.array([[y.min(), y.max()] for y,_ in axids]) # range accounting for panels
    xref = xrange[ref-1,:] # range for reference axes
    yref = yrange[ref-1,:]

    #-------------------------------------------------------------------------#
    # Get basemap.Basemap or cartopy.crs.Projection instances for map, and
    # override aspect ratio
    #-------------------------------------------------------------------------#
    # NOTE: Cannot have mutable dict as default arg, because it changes the
    # "default" if user calls function more than once! Swap dicts for None.
    proj = _notNone(projection, proj, None, names=('projection', 'proj'))
    proj_kw = _notNone(projection_kw, proj_kw, {}, names=('projection_kw', 'proj_kw'))
    proj    = _axes_dict(naxs, proj, kw=False, default='cartesian')
    proj_kw = _axes_dict(naxs, proj_kw, kw=True)
    basemap = _axes_dict(naxs, basemap, kw=False, default=False)
    axes_kw = {num:{} for num in range(1, naxs+1)}  # stores add_subplot arguments
    for num,name in proj.items():
        # The default, my CartesianAxes projection
        if name is None or name == 'cartesian':
            axes_kw[num]['projection'] = 'cartesian'
        # Builtin matplotlib polar axes, just use my overridden version
        elif name == 'polar':
            axes_kw[num]['projection'] = 'polar2'
            if num == ref:
                aspect = 1
        # Custom Basemap and Cartopy axes
        else:
            package = 'basemap' if basemap[num] else 'cartopy'
            obj, iaspect = projs.Proj(name, basemap=basemap[num], **proj_kw[num])
            if num == ref:
                aspect = iaspect
            axes_kw[num].update({'projection':package, 'map_projection':obj})

    #-------------------------------------------------------------------------#
    # Figure and axes panels
    #-------------------------------------------------------------------------#
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    # Create dictionaries of panel toggles and settings
    axpanels    = _notNone(axpanel, axpanels, '', names=('axpanel','axpanels'))
    axpanels_kw = _notNone(axpanel_kw, axpanels_kw, {}, names=('axpanel_kw', 'axpanels_kw'))
    axpanels    = _axes_dict(naxs, axpanels, kw=False, default='')
    axpanels_kw = _axes_dict(naxs, axpanels_kw, kw=True)

    # Get default figure panel keyword args
    # Also fill in empty keyword args
    kwout, kworig = {}, {}
    panel, panels = _notNone(panel, ''), _notNone(panels, '')
    allsides = panel + panels
    offsides = ''.join({*'lrbt'} - {*allsides})
    if len(allsides) != len({*allsides}):
        raise ValueError('You specified duplicate panel locations, e.g. panel="r" and colorbar="r".')
    keys = {*()}
    # Sanitize keyword args
    names = ('array', 'stack', 'share', 'space', 'width', 'sep')
    reg = re.compile(f'^[{offsides}]({"|".join(names)})$')
    ikw = {key:value for key,value in kwargs.items() if not reg.match(key)}
    ikwout, ikworig = _panels_kwargs(allsides, ikw, figure=True)
    keys.update(ikw.keys())
    kwout.update(ikwout)
    kworig.update(ikworig)
    # Warning
    jkw = {key:value for key,value in kwargs.items() if key not in keys}
    if jkw:
        warnings.warn(f'Ignoring figpanel keyword arg(s): {jkw}. Active sides are: {allsides!r}.')
    # Fill keyword args dicts -- unlike axes panels, entires are required
    for s in {*'lrbt'} - {*allsides}:
        kwout[s + 'array'] = None
        kwout[s + 'stack'] = 1
        kwout[s + 'share'] = False
        kwout[s + 'width'] = np.array([0])
        kwout[s + 'space'] = 0
        kwout[s + 'sep'] = np.array([])
        kworig[s + 'width'] = None
        kworig[s + 'space'] = None
        kworig[s + 'sep'] = None

    # Get figure panel arrays, depends on whether use passed 'panel' or 'panels'
    # which is not handled by _panels_kwargs!
    for ipanel,ispan in zip((panel,panels),(1,0)):
        for s in ipanel:
            nmax = (ncols if s in 'bt' else nrows)
            value = kwout[s + 'array']
            if value is None:
                value = ([1]*nmax if ispan else [*range(1,nmax+1)])
            elif not np.iterable(value) or len(value) != nmax:
                name = ('columns' if s in 'bt' else 'rows')
                raise ValueError(f'Need {nmax}-length list of integers for "{s}array" figure with {nmax} {name}, got {s}array={value!r}.')
            kwout[s + 'array'] = value

    #-------------------------------------------------------------------------#
    # Figure architecture
    #-------------------------------------------------------------------------#
    # Figure and/or axes dimensions
    names, values = (), ()
    if journal:
        figsize = _journals(journal) # if user passed width=<string > , will use that journal size
        spec = f'journal={journal!r}'
        names = ('axwidth', 'axheight', 'width')
        values = (axwidth, axheight, width)
        width, height = figsize
    elif figsize:
        spec = f'figsize={figsize!r}'
        names = ('axwidth', 'axheight', 'width', 'height')
        values = (axwidth, axheight, width, height)
        width, height = figsize
    elif width is not None or height is not None:
        spec = []
        if width is not None:
            spec.append(f'width={width!r}')
        if height is not None:
            spec.append(f'height={height!r}')
        spec = ', '.join(spec)
        names = ('axwidth', 'axheight')
        values = (axwidth, axheight)
    # Raise warning
    for name,value in zip(names,values):
        if value is not None:
            warnings.warn(f'You specified both {spec} and {name}={value!r}. Ignoring {name!r}.')
    # Standardize input
    width, height = units(width), units(height)
    axwidth, axheight = units(axwidth), units(axheight)

    # Store and standardize user input spacing args, so can override tight
    # layout calculations with them!
    left, right = units(left), units(right)
    bottom, top = units(bottom), units(top)
    # Space array with Nones allowed
    # Also store panel stack separations
    axwspace = np.atleast_1d(units(_notNone(wspace, space)))
    axhspace = np.atleast_1d(units(_notNone(hspace, space)))
    if len(axwspace) == 1:
        axwspace = np.repeat(axwspace, (ncols-1,))
    if len(axwspace) != ncols-1:
        raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(axwspace)}.')
    if len(axhspace) == 1:
        axhspace = np.repeat(axhspace, (nrows-1,))
    if len(axhspace) != nrows-1:
        raise ValueError(f'Require {nrows-1} height spacings for {nrows} rows, got {len(axhspace)}.')
    # Pull out originals and put them in arrays that match geometry
    # TODO: Why does wspace/hspace have to be an array?
    wspace = np.repeat((None,), 3*ncols + 1)
    hspace = np.repeat((None,), 3*nrows + 1)
    wspace[3:-1:3], hspace[3:-1:3] = axwspace, axhspace # between subplots, and leave panel slots empty so far
    wspace[0], wspace[-1] = kworig['lspace'], kworig['rspace'] # figure panels
    hspace[0], hspace[-1] = kworig['tspace'], kworig['bspace']
    wseps = np.repeat((None,), 3*ncols + 2)
    hseps = np.repeat((None,), 3*nrows + 2)
    wseps[0], wseps[-1] = kworig['lsep'], kworig['rsep']
    hseps[0], hseps[-1] = kworig['tsep'], kworig['bsep']
    wwidths = np.repeat((None,), 3*ncols + 2)
    hwidths = np.repeat((None,), 3*nrows + 2)
    wwidths[0], wwidths[-1] = kworig['lwidth'], kworig['rwidth']
    hwidths[0], hwidths[-1] = kworig['twidth'], kworig['bwidth']
    # Fill the subplots_orig_kw dictionary
    # NOTE: We add "seps" and "widths" *only* to subplots_orig_kw, to allow user
    # overrides -- info on current widths and seps is stored in ratios lists
    subplots_orig_kw = {
        'left':left, 'right':right, 'top':top, 'bottom':bottom,
        'wspace':wspace, 'hspace':hspace,
        'wwidths':wwidths, 'hwidths':hwidths, 'wseps':wseps, 'hseps':hseps,
        }

    # Default border spaces
    left   = _notNone(left,   units(rc['subplots.ylabspace']))
    right  = _notNone(right,  units(rc['subplots.innerspace']))
    top    = _notNone(top,    units(rc['subplots.titlespace']))
    bottom = _notNone(bottom, units(rc['subplots.xlabspace']))
    # Default spaces including panels
    axwspace[axwspace==None] = (
        units(rc['subplots.innerspace']) if sharey == 3
        else units(rc['subplots.ylabspace']) - units(rc['subplots.titlespace']) if sharey in (1,2) # space for tick labels only
        else units(rc['subplots.ylabspace']))
    axhspace[axhspace==None] = (
        units(rc['subplots.titlespace']) + units(rc['subplots.innerspace']) if sharex == 3
        else units(rc['subplots.xlabspace']) if sharex in (1,2) # space for tick labels and title
        else units(rc['subplots.titlespace']) + units(rc['subplots.xlabspace'])
        )
    wspace = np.repeat((None,), 3*ncols + 1)
    hspace = np.repeat((None,), 3*nrows + 1)
    wspace[3:-1:3], hspace[3:-1:3] = axwspace, axhspace # between subplots
    wspace[2:-1:3], hspace[2:-1:3] = 0, 0 # zero initial space for axes panels
    wspace[1:-1:3], hspace[1:-1:3] = 0, 0
    wspace[0], wspace[-1] = kwout['lspace'], kwout['rspace'] # figure panels
    hspace[0], hspace[-1] = kwout['tspace'], kwout['bspace']
    # Default ratios including panels
    axwidths = np.atleast_1d(_notNone(width_ratios, wratios, 1,
        names=('width_ratios', 'wratios')))
    axheights = np.atleast_1d(_notNone(height_ratios, hratios, 1,
        names=('height_ratios', 'hratios')))
    if len(axwidths) == 1:
        axwidths = np.repeat(axwidths, (ncols,))
    if len(axheights) == 1:
        axheights = np.repeat(axheights, (nrows,))
    if len(axwidths) != ncols:
        raise ValueError(f'Got {ncols} columns, but {len(axwidths)} wratios.')
    if len(axheights) != nrows:
        raise ValueError(f'Got {nrows} rows, but {len(axheights)} hratios.')
    wratios = np.repeat((None,), 3*ncols + 2)
    hratios = np.repeat((None,), 3*nrows + 2)
    wratios[2:-1:3], hratios[2:-1:3] = axwidths, axheights # subplots
    wratios[3:-1:3], hratios[3:-1:3] = 0, 0 # zero width initial axes panels
    wratios[1:-1:3], hratios[1:-1:3] = 0, 0
    wratios[0]  = sum(kwout['lwidth']) + sum(kwout['lsep'])
    wratios[-1] = sum(kwout['rwidth']) + sum(kwout['lsep'])
    hratios[0]  = sum(kwout['twidth']) + sum(kwout['tsep'])
    hratios[-1] = sum(kwout['bwidth']) + sum(kwout['bsep'])

    # Parse arguments, fix dimensions in light of desired aspect ratio
    figsize, gridspec_kw, subplots_kw = _subplots_geometry(
        nrows=nrows, ncols=ncols,
        aspect=aspect, xref=xref, yref=yref,
        left=left, right=right, bottom=bottom, top=top,
        width=width, height=height, axwidth=axwidth, axheight=axheight,
        wratios=wratios, hratios=hratios, wspace=wspace, hspace=hspace,
        )
    fig = plt.figure(FigureClass=Figure, tight=tight, figsize=figsize, ref=ref,
        pad=pad, axpad=axpad, panelpad=panelpad, autoformat=autoformat,
        subplots_orig_kw=subplots_orig_kw, subplots_kw=subplots_kw,
        gridspec_kw=gridspec_kw)
    gridspec = fig._main_gridspec

    #-------------------------------------------------------------------------#
    # Draw on figure
    #-------------------------------------------------------------------------#
    # Draw main subplots
    axs = naxs*[None] # list of axes
    for idx in range(naxs):
        # Get figure gridspec ranges
        num = idx + 1
        x0, x1 = xrange[idx,0]*3 + 2, xrange[idx,1]*3 + 2 # offsets for panels
        y0, y1 = yrange[idx,0]*3 + 2, yrange[idx,1]*3 + 2
        # Draw subplot
        subplotspec = gridspec[y0:y1+1, x0:x1+1]
        with fig._unlock():
            axs[idx] = fig.add_subplot(subplotspec, number=num,
                spanx=spanx, spany=spany, alignx=alignx, aligny=aligny,
                sharex_level=sharex, sharey_level=sharey,
                **axes_kw[num])

    # Set up shared axes and assign main axes
    fig._main_axes = axs
    fig._share_axes_setup()

    # Draw axes panels after all main subplots are drawn
    # NOTE: This *must* come after shared axes are set up! Otherwise tight
    # layout scaling is wrong
    names = ('stack', 'share', 'space', 'width', 'sep')
    for idx in range(naxs):
        ax = axs[idx]
        num = idx + 1
        keys = {*()}
        sides, panels_kw = axpanels[num], axpanels_kw[num]
        # Loop through sides
        for s in sides:
            offsides = ''.join({*'lrbt'} - {s})
            reg = re.compile(f'^[{offsides}]({"|".join(names)})$')
            ikw = {key:value for key,value in panels_kw.items() if not reg.match(key)}
            keys.update(ikw.keys())
            fig._add_panel(ax, s, order=order, **ikw)
        # Warning message
        jkw = {key:value for key,value in panels_kw.items() if key not in keys}
        if jkw:
            warnings.warn(f'Ignoring axpanels_kw keyword arg(s): {jkw}. Active sides are: axpanels={sides!r}.')

    # Draw figure panels
    for s in 'blrt':
        # Get keyword args and gridspec args
        array = kwout[s + 'array']
        width = kwout[s + 'width']
        share = kwout[s + 'share']
        sep = kwout[s + 'sep']
        if array is None:
            continue
        if s in 'lr':
            hspace, hratios, wspace, wratios = [], [1], sep, width
        else:
            hspace, hratios, wspace, wratios = sep, width, [], [1]
        nrows, ncols = len(hratios), len(wratios)
        # Loop through unique numbers
        paxs = []
        side = _side_translate[s]
        for num in np.unique(array).flat:
            # Get subplotspec and insert gridspec
            if num == 0:
                continue
            idx, = np.where(array == num)
            idx1 = slice(min(idx)*3 + 2, max(idx)*3 + 3) # ignores axes panel slots
            idx2 = 0 if s in 'lt' else -1
            if s in 'bt':
                idx1, idx2 = idx2, idx1
            subplotspec = gridspec[idx1,idx2]
            igridspec = FlexibleGridSpecFromSubplotSpec(
                subplot_spec=subplotspec,
                nrows=nrows, ncols=ncols,
                wspace=wspace, hspace=hspace,
                width_ratios=wratios, height_ratios=hratios,
                )
            # Draw axes
            ipaxs = []
            for i in range(max((nrows,ncols))):
                with fig._unlock():
                    ipax = fig.add_subplot(igridspec[i],
                        projection='panel',
                        side=side, share=share,
                        )
                ipaxs += [ipax]
            paxs += [ipaxs]
        # Sort panel axes into row-major or column-major order
        if (s in 'bt' and order == 'C') or (s in 'lr' and order != 'C'):
            paxs = [*zip(*paxs)]
        # Store in axes_grid with support for 2D indexing
        n = len(paxs[0])
        paxs = [ax for ipaxs in paxs for ax in ipaxs] # unfurl
        paxs = axes_grid(paxs, n=n, order=order)
        setattr(fig, '_' + side + 'panel', paxs)

    # Return figure and axes
    n = (ncols if order == 'C' else nrows)
    return fig, axes_grid(axs, n=n, order=order)

