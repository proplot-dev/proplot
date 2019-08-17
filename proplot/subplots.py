#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures and axes.
The `subplots` function is all you'll need to directly use here.
It returns a `Figure` instance and an `axes_grid` container of
`~proplot.axes.BaseAxes` axes, whose positions are controlled by the
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
import os
import re
import numpy as np
import functools
import warnings
import matplotlib.pyplot as plt
import matplotlib.figure as mfigure
import matplotlib.transforms as mtransforms
import matplotlib.gridspec as mgridspec
from .rctools import rc
from .utils import _default, _counter, _timer, units
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
_side_translate_inverse = {
    'left':   'l',
    'right':  'r',
    'bottom': 'b',
    'top':    't',
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
        """Wraps the string representation."""
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
        .. code-block:: python

            import proplot as plot
            f, axs = plot.subplots(nrows=3, ncols=3, colorbars='b', bstack=2)
            axs[0] # the subplot in the top-right corner
            axs[3] # the first subplot in the second row
            axs[1,2] # the subplot in the second row, third from the left
            axs[:,0] # the subplots in the first column

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
        .. code-block:: python

            import proplot as plot
            f, axs = plot.subplots(nrows=2, ncols=2, axcolorbars='b')
            axs.format(xtick=5) # calls "format" on all subplots in the list
            axs.bpanel.colorbar(m) # calls "colorbar" on all panels in the axes_grid returned by "axs.bpanel"

        """
        attrs = (*(getattr(ax, attr, None) for ax in self),)
        # Not found
        if None in attrs:
            raise AttributeError(f'Attribute "{attr}" not found.')
        # Empty
        if not attrs:
            def null_iterator(*args, **kwargs):
                return None
            return null_iterator
        # Panels
        if all(isinstance(_, (axes_grid, axes.BaseAxes, axes.EmptyPanel))
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
        raise AttributeError(f'Found mixed types for attribute "{attr}".')

#-----------------------------------------------------------------------------#
# Gridspec classes
#-----------------------------------------------------------------------------#
def _adjust(n):
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
            subplots. Values are scaled relative to the height and width
            ratios. For example, ``wspace=0.1`` yields a space 10% the width
            of the axes.

            If list, length of ``wspace``, ``hspace`` must be ``ncols-1``,
            ``nrows-1``.
        height_ratios, width_ratios : list of float
            Ratios for the width/height of columns/rows of subplots.
            For example, ``width_ratios=[1,2]`` specifes 2 columns of subplots,
            the second one twice as wide as the first.
        left, right, top, bottom : float or str
            Passed to `~matplotlib.gridspec.GridSpec`. Indicates width of "margins"
            surrounding the grid. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.

            Generally, these are used to set the "figure edges" around the
            region of subplots. If `~proplot.subplots.subplots` was called
            with ``tight=True`` (the default), these are ignored.
        """
        # Add these as attributes; want _spaces_as_ratios to be
        # self-contained, so it can be invoked on already instantiated
        # gridspec (see update() method).
        # TODO Does _nrows or _ncols conflict with default gridspec attributes?
        self._nrows = nrows*2-1 # used for get_geometry and needed by _spaces_as_ratios
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
                raise ValueError(f'Invalid index: "{key}".')
            num1 = _normalize(k1, nrows_visible)
            num2 = _normalize(k2, ncols_visible)
            num1, num2 = np.ravel_multi_index((num1, num2), (nrows, ncols))
        # Correct for negative nums
        num1, num2 = _adjust(num1), _adjust(num2)
        return mgridspec.SubplotSpec(self, num1, num2)

    def _spaces_as_ratios(self,
            hspace=None, wspace=None, # spacing between axes
            height_ratios=None, width_ratios=None,
            **kwargs):
        """For keyword arg usage, see `FlexibleGridSpecBase`."""
        # Parse flexible input
        nrows, ncols = self.get_visible_geometry()
        hratios = np.atleast_1d(_default(height_ratios, 1))
        wratios = np.atleast_1d(_default(width_ratios,  1))
        hspace = np.atleast_1d(_default(hspace, np.mean(hratios)*0.10)) # this is relative to axes
        wspace = np.atleast_1d(_default(wspace, np.mean(wratios)*0.10))
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
        if ncols > 1 and len(wspace) != ncols-1:
            raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(wspace)}.')
        if nrows > 1 and len(hspace) != nrows-1:
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

    def update(self, **kwargs):
        """Updates the width ratios, height ratios, and spacing for subplot
        columns and rows."""
        wratios, hratios, kwargs = self._spaces_as_ratios(**kwargs)
        self.set_width_ratios(wratios)
        self.set_height_ratios(hratios)
        for key in ('nrows','ncols'):
            kwargs.pop(key, None) # cannot be modified
        super().update(**kwargs) # remaining kwargs should just be left, right, top, bottom

    def get_visible_geometry(self):
        """Like `~matplotlib.gridspec.GridspecBase.get_geometry`, but returns
        the number of visible rows and columns, i.e. the number of rows and
        columns that aren't skipped over by `~FlexibleGridSpecBase.__getitem__`."""
        return self._nrows_visible, self._ncols_visible

class FlexibleGridSpec(FlexibleGridSpecBase, mgridspec.GridSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpec`."""
    pass

class FlexibleGridSpecFromSubplotSpec(FlexibleGridSpecBase, mgridspec.GridSpecFromSubplotSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpecFromSubplotSpec`."""
    pass

#-----------------------------------------------------------------------------#
# Figure class
#-----------------------------------------------------------------------------#
def _xtight(ax):
    """Returns tight bbox coordinates."""
    bbox = ax._tight_bbox
    return (bbox.xmin, bbox.xmax)

def _ytight(ax):
    """Returns tight bbox coordinates."""
    bbox = ax._tight_bbox
    return (bbox.ymin, bbox.ymax)

def _xrange(ax):
    """Gets the column range for the axes. If this is a child of the main axes,
    returns properties for the parent axes."""
    # Factor of two corrects for inaccessible 'space' gridspec locations
    # e.g. 0-- > 0, 2-- > 1, 4-- > 2, and axes can never lie on 1, 3, 5, etc.
    subplotspec = ax.get_subplotspec().get_topmost_subplotspec()
    _, _, _, _, col1, col2 = subplotspec.get_rows_columns()
    return col1//2, col2//2

def _yrange(ax):
    """Gets the row range for axes. If this is a child of the main axes,
    returns properties for the parent axes."""
    subplotspec = ax.get_subplotspec().get_topmost_subplotspec()
    _, _, row1, row2, _, _ = subplotspec.get_rows_columns()
    return row1//2, row2//2

def _xtight_full(ax, children=True):
    """Gets *x* span of tight bounding box, including panels, and shared axes."""
    # TODO: Better resting for axes visibility
    axs = ax._iter_children_panels()
    axs = [ax for ax in axs if ax._tight_bbox is not None]
    xs = np.array([_xtight(ax) for ax in axs])
    return xs[:,0].min(), xs[:,1].max()

def _ytight_full(ax, children=True):
    """Get span, accounting for panels, shared axes, and whether axes has
    been replaced by colorbar in same location."""
    axs = ax._iter_children_panels()
    axs = [ax for ax in axs if ax._tight_bbox is not None]
    ys = np.array([_ytight(ax) for ax in axs])
    return ys[:,0].min(), ys[:,1].max()

class Figure(mfigure.Figure):
    def __init__(self,
        tight=None, tightborders=None, tightsubplots=None, tightpanels=None,
        flush=False, wflush=None, hflush=None,
        borderpad=None, subplotpad=None, panelpad=None,
        autoformat=True, ref=1, # ref is documented in subplots
        main_gridspec=None, subplots_kw=None,
        **kwargs):
        """
        The `~matplotlib.figure.Figure` instance returned by `subplots`. At
        draw-time, an improved "tight layout" adjustment is applied, and
        the space around the figure edge, between subplots, and between
        panels is changed to accommodate subplot content. Figure dimensions
        may be automatically scaled to preserve subplot aspect ratios.

        Parameters
        ----------
        tight : bool, optional
            If passed, sets `tightborders`, `tightsubplots`, and `tightpanels`
            all at once.
        tightborders : bool, optional
            Whether to draw a tight bounding box around the whole figure.
            Defaults to `tight` if passed, ``rc['tight']`` otherwise.
        tightsubplots : bool, optional
            Whether to automatically space out subplots to prevent overlapping
            axis tick labels, etc. Defaults to `tight` if passed,
            ``rc['tight']`` otherwise.
        tightpanels : bool, optional
            Whether to automatically space between subplots and their panels
            to prevent overlap.  Defaults to `tight` if passed,
            ``rc['tight']`` otherwise.
        borderpad : float or str, optional
            Margin size for tight bounding box surrounding the edge of the
            figure. Defaults to ``rc['subplots.borderpad']``. If float, units are inches.
            If string, units are interpreted by `~proplot.utils.units`.
        subplotpad : float or str, optional
            Margin size between content from adjacent subplots. Defaults to
            ``rc['subplots.subplotpad']``. If float, units are inches.
            If string, units are interpreted by `~proplot.utils.units`.
        panelpad : float or str, optional
            Margin size between content from subplots and their child panels.
            Defaults to ``rc['subplots.subplotpad']``. If float, units are inches.
            If string, units are interpreted by `~proplot.utils.units`.
        flush, wflush, hflush : bool, optional
            Whether subplots should be "flush" against each other in the
            horizontal (`wflush`), vertical (`hflush`), or both (`flush`)
            directions. Useful if you want to let axis ticks overlap with other
            axes, and just want axes spines touching each other. Note that
            instead of ``flush=0``, you can also use ``tightsubplots=False``
            with manual gridspec spacings ``wspace=0`` and ``hspace=0``.
        autoformat : bool, optional
            Whether to automatically format the axes when a `~pandas.Series`,
            `~pandas.DataFrame` or `~xarray.DataArray` is passed to a plotting
            command.
        main_gridspec, subplots_kw
            The primary gridspec used for populating the figure with
            subplots, and the keyword args passed to the `subplots` command
            and used to generate this figure.

        Other parameters
        ----------------
        **kwargs
            Passed to `matplotlib.figure.Figure`.

        Warning
        -------
        To make "tight layout" and aspect ratio adjustments, we rely on the
        hidden `~Figure._update_layout` method
        `~Figure._update_layout` may also call
        `~matplotlib.figure.Figure.set_size_inches`, which may interfere with
        figure rendering when used with certain popup backends. It also seems
        to have
        Inline
        backends seem to have no problem with this.
        """
        # Tight settings
        # NOTE: _tight taken by matplotlib tight layout, use _smart_tight!
        self._tight_borders = _default(tightborders, tight, rc['tight'])
        self._tight_subplots = _default(tightsubplots, tight, rc['tight'])
        self._tight_panels = _default(tightpanels, tight, rc['tight'])
        self._smart_tight = (self._tight_borders or self._tight_subplots or self._tight_panels)
        self._post_init = True
        self._smart_tight_init = True
        self._border_pad = units(_default(borderpad, rc['subplots.borderpad']))
        self._subplot_pad  = units(_default(subplotpad,  rc['subplots.subplotpad']))
        self._panel_pad = units(_default(panelpad, rc['subplots.panelpad']))
        self._wflush = _default(wflush, flush)
        self._hflush = _default(hflush, flush)
        self._locked = True
        self._autoformat = autoformat
        # Gridspec info and panels, filled in by subplots
        # TODO: Use existing axes tracking attributes?
        self._ref_num = ref
        self._main_axes = []
        self._main_gridspec = main_gridspec
        self._axes_spanning = []
        self._subplots_kw = subplots_kw
        self._leftpanel   = axes.EmptyPanel()
        self._bottompanel = axes.EmptyPanel()
        self._rightpanel  = axes.EmptyPanel()
        self._toppanel    = axes.EmptyPanel()
        # Initialize
        super().__init__(**kwargs)
        self.suptitle('') # add _suptitle attribute

    def _get_tight_bboxs(self, renderer=None):
        """Sets the ``_tight_bbox`` attribute on axes, so that we don't have
        to call the tight bbox algorithm multiple times."""
        for ax in (*self._main_axes, *self.leftpanel, *self.bottompanel,
                   *self.rightpanel, *self.toppanel):
            if not ax:
                continue
            for iax in ax._iter_children_panels():
                bbox = iax.get_tightbbox(renderer)
                iax._tight_bbox = bbox

    def _get_side_axes(self, side, ref=None):
        """Returns groups of axes in row or column or the single group in the
        same row or column as axes `ref`."""
        idx = (1 if side in 'br' else 0) # which side of range to test
        _along = (_yrange if side in 'bt' else _xrange)
        if ref is None:
            nrows, ncols = self._main_gridspec.get_visible_geometry()
            nums = [*range(nrows if side in 'bt' else ncols)]
        else:
            ref = _along(ref)[idx] # side for a particular axes
            nums = [ref]
        groups = []
        for num in nums:
            axs = [ax for ax in self._main_axes if _along(ax)[idx] == num]
            groups.append(axs)
        if ref is None:
            return groups
        else:
            return groups[0]

    def _get_shared_axes(self, mode, ref=None):
        """Returns groups of axes that share the same horizontal or vertical
        extent, used for determining shared axes. Input `mode` should be
        ``'x'`` or ``'y'``. Can return list of groups or just the single
        group matching subplot `ref`."""
        # This does not already exist in matplotlib API. Matplotlib applies
        # automatic sharing for simple grids, not complex ones.
        idx = (1 if mode == 'x' else 0)
        base_idx = (np.argmax if mode == 'x' else np.argmin)
        _along = (_xrange if mode == 'x' else _yrange)
        _across = (_yrange if mode == 'x' else _xrange) # for x-sharing, want axes that span same extent in y direction
        if ref is None:
            axs = self._main_axes
            ranges = {tuple(_along(ax)) for ax in axs} # unique ranges
        else:
            axs = [ref]
            ranges = {tuple(_along(ref))} # unique ranges, just the one
        groups = []
        for irange in ranges:
            iaxs = [ax for ax in axs if tuple(_along(ax)) == irange]
            ibase = iaxs.pop(base_idx([_across(iax)[idx] for iax in iaxs]))
            groups.append([ibase, *iaxs])
        if ref is None:
            return groups
        else:
            return groups[0]

    def _get_space_axes(self, mode, ref=None):
        """Returns groups of axes that abutt against the same space in a given
        column or row, used for tight layout adjustment. Input `mode` should
        be ``'x'`` or ``'y'``. Output is list of lists of length-2 lists, where
        the length-2 lists in each sublist are unique groups of axes on that
        particular interface whose right, left (``mode='x'``) or top, bottom
        (``mode='y'``) *edges* abutt against each other. Breaking down into
        these groups keeps us from applying unnecessary padding if e.g. the
        axes in row 1, column 1 has a large xlabel and the axes in row 2,
        column 2 has a large title."""
        # Initial stuff
        nrows, ncols = self._main_gridspec.get_visible_geometry()
        nspaces = (ncols if mode == 'x' else nrows)
        nacross = (nrows if mode == 'x' else ncols)
        _along = (_xrange if mode == 'x' else _yrange)
        _across = (_yrange if mode == 'x' else _xrange)
        if ref is not None:
            raise NotImplementedError
        axs = [ax for ax in self._main_axes if ax.get_visible()]
        if self.leftpanel:
            axs.extend(self.leftpanel[:,-1])
        if self.rightpanel:
            axs.extend(self.rightpanel[:,0])
        if self.bottompanel:
            axs.extend(self.bottompanel[0,:])
        if self.toppanel:
            axs.extend(self.toppanel[-1,:])
        # Iterate along spaces
        groups = []
        ranges_along = np.array([_along(ax) for ax in axs])
        ranges_across = np.array([_across(ax) for ax in axs])
        for space in range(1, nspaces): # spaces between e.g. columns
            igroups = []
            filt1 = ranges_along[:,1]+1 == space # i.e. right/bottom edge abutts against this space
            filt2 = ranges_along[:,0] == space # i.e. left/top edge abutts against this space
            for idx in range(nacross): # e.g. each row
                # Find axes that abutt aginst this space on this particular row
                filt = (ranges_across[:,0] <= idx) & (idx <= ranges_across[:,1]) # i.e. subplot spans this particular column or row
                if sum(filt) < 2: # no interface here
                    continue
                idx1, = np.where(filt & filt1)
                idx2, = np.where(filt & filt2)
                if not (idx1.size == 1 and idx2.size == 1):
                    continue # normally both zero, but one can be non-zero e.g. if have left and bottom panel, empty space in corner
                # Put these axes into unique groups. Store groups as
                # (left axes, right axes) or (bottom axes, top axes) pairs.
                new = True
                ax1, ax2 = axs[idx1[0]], axs[idx2[0]]
                if mode != 'x':
                    ax1, ax2 = ax2, ax1 # yrange is top-to-bottom, so make this bottom-to-top
                for group in igroups:
                    if ax1 in group[0] or ax2 in group[1]:
                        group[0].update({ax1})
                        group[1].update({ax2})
                        new = False
                        break
                if new:
                    igroups.append([{ax1}, {ax2}]) # form new group
            groups.append(igroups)
        return groups

    @_timer
    def _add_panel(self, ax, side, order='C', mode='panel', **kwargs):
        """Hidden method that powers `~proplot.axes.panel_axes`. Makes more sense
        to define in subplots, because it does a bunch of alignment steps
        that rely on the figure instance."""
        # Checks
        side = _side_translate_inverse.get(side, side) # single char
        if side not in (*'lrbt',):
            raise ValueError(f'Invalid side "{side}".')
        name = _side_translate[side] # full word
        paxs = getattr(ax, name + 'panel')
        if paxs:
            warnings.warn(f'{name}panel already exists.')
            return paxs # return existing panels! should already be synced
        gridspec = ax._main_gridspec
        if gridspec is None:
            raise ValueError(f'Gridspec attribute is empty.')

        # Parse keyword args
        # NOTE: This looks funny but it is meant to work with alternative API
        # where axes panels are added all at once!
        if mode=='panel':
            args, kwtuple = (side, '', ''), (kwargs, {}, {})
        elif mode=='colorbar':
            args, kwtuple = ('', side, ''), ({}, kwargs, {})
        elif mode=='legend':
            args, kwtuple = ('', '', side), ({}, {}, kwargs)
        else:
            raise ValueError(f'Invalid panel mode "{mode}". Must be one of "panel", "colorbar", or "legend".')
        kwargs = _panels_kwargs(*args, *kwtuple, figure=False)
        kwargs.pop('sides') # which is only needed for axpanels='' usage
        # Get props pertaining to this side, issue warning if found others!
        width, sep = kwargs.pop(side + 'width'), kwargs.pop(side + 'sep')
        space, flush = kwargs.pop(side + 'space'), kwargs.pop(side + 'flush')
        share = kwargs.pop(side + 'share')
        for key in (*kwargs,): # so dict can change size
            if kwargs[key] is None:
                kwargs.pop(key)
        if kwargs:
            raise ValueError(f'Unknown or unused keyword args for {name}panel: {kwargs}.')

        # Fix main gridspec
        if side in 'lr':
            ratios = gridspec.get_width_ratios()
        else:
            ratios = gridspec.get_height_ratios()
        idx = (1 if side in 'lt' else -2)
        ratios[idx] = space # space between subplot and panels
        idx = (0 if side in 'lt' else -1)
        ratios[idx] = sum(sep) + sum(width)
        # Fix axes ratio, so it is in physical units
        # ratios[2] -= (sum(sep) + sum(width) + space)
        # Get subplotspec
        if side == 'l':
            key = (1,0)
        elif side == 'r':
            key = (1,2)
        elif side == 'b':
            key = (2,1)
        elif side == 't':
            key = (0,1)
        subplotspec = gridspec[key]
        # Stack settings
        stack = len(width)
        if flush:
            sep = sep*0.0
        if side in 'lr':
            wspace, hspace = sep, []
            wratios, hratios = width, [1]
            nrows, ncols = 1, stack
        else:
            wspace, hspace = [], sep
            wratios, hratios = [1], width
            nrows, ncols = stack, 1

        # Draw panel(s)
        paxs = []
        sgridspec = FlexibleGridSpecFromSubplotSpec(
            subplot_spec=subplotspec,
            nrows=nrows, ncols=ncols,
            wspace=wspace, hspace=hspace,
            width_ratios=wratios, height_ratios=hratios)
        for i in range(stack):
            self._locked = False
            pax = self.add_subplot(sgridspec[i],
                side=name, flush=flush,
                share=share, parent=ax, projection='panel',
                )
            self._locked = True
            pax._main_gridspec = gridspec
            pax._stack_gridspec = sgridspec
            paxs += [pax]

        # Add as axes_grid. Support 2D indexing, even though these are
        # always vector stacks, because consistency. See axes_grid docs.
        n = 1 if (side in 'tb' and order == 'C') or (side in 'lr' and order != 'C') else stack
        paxs = axes_grid(paxs, n=n, order=order)
        setattr(ax, '_' + name + 'panel', paxs) # hidden attribute for property
        self._sync_panels(ax, side, sep=sep, width=width, space=space)
        return paxs

    def _align_suplabels(self, renderer):
        """Adjusts position of row titles and figure super title."""
        # Adjust row labels. Remember we always leave room for figure panels
        # gridspec, so want to find where xrange and yrange == 1, not 0.
        # NOTE: Must uset get_tightbbox so this will still work if we are
        # not doing tight layout mode, and because bbox coordinates may be
        # invalid after figure has been resized in _smart_tight_layout.
        width, height = self.get_size_inches()
        axs = [ax for ax in self._main_axes if _xrange(ax)[0] == 1] # order by xrange
        lxs, labels = [], []
        for ax in axs:
            label = ax.rowlabel
            if not label.get_text().strip():
                continue
            # Iterate panels
            ixs = []
            label.set_visible(False) # make temporarily invisible, so tightbbox does not include existing label!
            if ax.leftpanel:
                iaxs = ax.leftpanel
            else:
                iaxs = (ax, *ax.toppanel, *ax.bottompanel)
            for iax in iaxs:
                if not iax:
                    continue
                for jax in iax._iter_children():
                    # Box for column label
                    bbox = jax.get_tightbbox(renderer)
                    x, _ = self.transFigure.inverted().transform((bbox.xmin, 0))
                    ixs.append(x)
            label.set_visible(True)
            # Update position
            if not ixs:
                warnings.warn('Axes on left row are invisible. Cannot determine rowtitle position.')
                continue
            x = min(ixs) - (0.6*label.get_fontsize()/72)/width # see: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
            lxs.append(x)
            labels.append(label)
        # Apply
        if lxs:
            x = min(lxs)
            for label in labels:
                label.update({'x':x, 'y':0.5})

        # Adjust col labels -- this is much simpler
        sys, saxes = [], [] # for suptitle
        lys, labels = [], []
        suptitle = self._suptitle
        suptitle.set_visible(False)
        axs = [ax for ax in self._main_axes if _yrange(ax)[0] == 1] # order by xrange
        for ax in axs:
            label = ax.collabel
            if not label.get_text().strip() and not suptitle.get_text().strip():
                continue
            # Get maximum y-position among all children, necessary because
            # e.g. twin x-axes with their own labels are common
            iys = []
            label.set_visible(False)
            if ax.toppanel and ax.toppanel.get_visible():
                iaxs = ax.toppanel
            else:
                iaxs = (ax, *ax.leftpanel, *ax.rightpanel)
            for iax in iaxs:
                if not iax:
                    continue
                for jax in iax._iter_children():
                    # Box for column label
                    bbox = jax.get_tightbbox(renderer)
                    _, y = self.transFigure.inverted().transform((0, bbox.ymax))
                    iys.append(y)
            label.set_visible(True)
            # Update column label position
            if not iys:
                warnings.warn('Axes on top row is invisible. Cannot determine coltitle position.')
                continue
            saxes.append(ax)
            if label.get_text().strip():
                y = max(iys) + (0.3*label.get_fontsize()/72)/height # see: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
                sys.append(None) # we have to re-calcualte
                lys.append(y)
                labels.append(label)
            else:
                sys.append(max(iys)) # just use the coordinate we were going to use for the label
        # Apply
        if lys:
            y = max(lys)
            for label in labels:
                label.update({'x':0.5, 'y':y})

        # Super title position
        # Get x position as center between left edge of leftmost main axes
        # and right edge of rightmost main axes
        suptitle.set_visible(True)
        if suptitle.get_text().strip():
            iys = []
            for y,ax in zip(sys,saxes):
                if y is None:
                    bbox = ax.get_tightbbox(renderer)
                    _, y = self.transFigure.inverted().transform((0, bbox.ymax))
                iys.append(y)
            if not iys: # means sys and saxes is empty
                warnings.warn('All axes on top row are invisible. Cannot determine suptitle position.')
            else:
                kw = self._subplots_kw
                lsep, rsep = _default(kw['lsep'], []), _default(kw['rsep'], [])
                lwidth = _default(kw['lwidth'], [0])
                rwidth = _default(kw['rwidth'], [0])
                lspace = _default(kw['lspace'], 0)
                rspace = _default(kw['rspace'], 0)
                left = kw['left'] + sum(lwidth) + sum(lsep) + lspace
                right = kw['right'] + sum(rwidth) + sum(rsep) + rspace
                x = left/width + 0.5*(width - left - right)/width
                y = max(iys) + (0.3*suptitle.get_fontsize()/72)/height
                suptitle.update({'x':x, 'y':y, 'ha':'center', 'va':'bottom', 'transform':self.transFigure})

    def _update_layout(self, renderer=None):
        """Aligns row labels, column labels, and super titles, and applies
        tight layout automatic spacing."""
        if renderer is None:
            renderer = self.canvas.get_renderer()
        # Just align subplot labels
        if not self._smart_tight:
            self._align_suplabels(renderer)
        # Align labels and apply tight spacing
        else:
            # Get boxes before aligning labels -- get_tightbbox will
            # call title offset post-processing step
            for ax in self._main_axes:
                ax.rowlabel.set_visible(False)
                ax.collabel.set_visible(False)
            self._get_tight_bboxs(renderer) # can use same bboxs
            for ax in self._main_axes:
                ax.rowlabel.set_visible(True)
                ax.collabel.set_visible(True)
            # Update gridspec props and label positions
            # Must align figlabels twice! Cannot just align subplots, then
            # align labels, then get tight borders. Dunno why it fails.
            self._align_suplabels(renderer)
            self._smart_tight_layout(renderer)
            self._align_suplabels(renderer)
        # Flag
        self._post_init = False

    def _panel_tight_layout(self, side, axs, renderer, figure=False):
        """From a list of axes in the same row or column, figure out the
        necessary new 'spacing' to apply to the panel axes gridspec object.
        For axes panels, this function just modifies mutable width and height
        ratio lists in place, and returns the 'wpanels' or 'hpanels' arguments.
        For figure panels, it returns the 'lsep', 'rsep', or 'bsep' argument
        needed for _subplots_geometry."""
        # Initial stuff
        # TODO: Make sure this works for all combinations of figure panels!
        # NOTE: We edit the *ratio lists* themselves. The getters do not
        # return copies, but the mutable underlying objects!
        pad = self._panel_pad
        if all(isinstance(ax, (axes_grid, axes.PanelAxes, axes.EmptyPanel)) for ax in axs):
            paxs = axs
            paxs_main = [pax[-1] if side in 'tl' else pax[0] for pax in paxs]
            gridspecs = [pax._main_gridspec if pax else None
                            for pax in paxs_main]
        else:
            paxs = [getattr(ax, side + 'panel') for ax in axs]
            paxs_main = [pax[-1] if side in 'tl' else pax[0] for pax in paxs]
            gridspecs = [ax._main_gridspec for ax in axs]
        if not any(paxs): # all are EmptyPanel
            return 0
        elif len({len(pax) for pax in paxs if pax}) > 1:
            warnings.warn(f'Panel tight layout failed on side {side}. Conflicting number of stacked panels in same row or column: {[len(pax) for pax in paxs]}.')
        # Outer gridspec ratios
        if side in 'lr':
            ratios = [gs.get_width_ratios() if gs else None
                        for gs in gridspecs]
        else:
            ratios = [gs.get_height_ratios() if gs else None
                        for gs in gridspecs]
        # Inner 'stack' gridspec ratios
        pgridspecs = [pax._stack_gridspec if pax else None
                        for pax in paxs_main]
        if side in 'lr':
            pratios = [gs.get_width_ratios() if gs else None
                        for gs in pgridspecs]
        else:
            pratios = [gs.get_height_ratios() if gs else None
                        for gs in pgridspecs]

        # Get necessary space between stacked panels, iterate through pairs
        # of panels and get *minimum* actual space. Space is calculated:
        # 1) Bottom of top panel minus top of bottom panel
        # 2) Left of right panel minus right of left panel
        # NOTE: Order of stacks is always left-right and top-bottom
        seps = []
        for ax,pax in zip(axs,paxs):
            isep = []
            for i in range(len(pax)-1): # empty if pax is EmptyPanel
                ipaxs = pax[i:i+2] # lists are left-to-right, top-to-bottom
                if not all(ipaxs): # found EmptyPanel
                    continue
                if side in 'lr':
                    ispans = [_xtight(pax) for pax in ipaxs]
                    isep.append((ispans[1][0] - ispans[0][1])/self.dpi) # bottom of top one minus top of bottom one
                else: # 'tb'
                    ispans = [_ytight(pax) for pax in ipaxs]
                    isep.append((ispans[0][0] - ispans[1][1])/self.dpi) # bottom of top one minus top of bottom one
            seps.append(isep)

        # Group seps so they correspond to like rows or columns
        # If seps is [[], ..., []], results will be []
        seps = [sep for sep in seps if sep] # ignore EmptyPanel
        seps = [[*_] for _ in zip(*seps)]
        # Apply updated seps to the underlying panel gridspec and main
        # axes gridspec width and height ratios
        sep = []
        flush = False
        if any(paxs_main):
            flush = [pax for pax in paxs_main if pax][0]._flush
        for idx,isep in enumerate(seps):
            idx = 1 + idx*2 # index of *spaces* in ratios list
            sep.append(
                0 if flush else max([0, pratios[0][idx] - min(isep) + pad])
                )
        for iratios,ipratios in zip(ratios,pratios):
            if not iratios and not figure:
                raise RuntimeError(f'This should not happen with axes panels.')
            if not ipratios:
                continue
            # Adjust inner gridspec ratios
            for idx,isep in enumerate(sep):
                idx = 1 + idx*2 # index of *spaces* in ratios list
                ipratios[idx] = isep
            # Adjust outer gridspec ratios
            if figure:
                continue
            idx_stack = (0 if side in 'lt' else -1)
            iratios[idx_stack] = sum(ipratios)

        # Conditional return
        if figure: # just need the 'sep' argument
            return sep
        else:
            # Get space between main subplot and adjacent panels. Iterate
            # through panels, get maximum necessary spacing.
            spaces = []
            for pax in paxs_main:
                if not pax: # ignore EmptyPanel
                    continue
                if side in 'lr':
                    pspan = _xtight(pax)
                    span = _xtight(pax._parent)
                else:
                    pspan = _ytight(pax)
                    span = _ytight(pax._parent)
                if side in 'tr':
                    spaces.append((pspan[0] - span[1])/self.dpi)
                else:
                    spaces.append((span[0] - pspan[1])/self.dpi)
            # Assign new spacing to ratios
            idx = (-2 if side in 'br' else 1)
            space = (0 if flush or not spaces else
                     max(0, ratios[0][idx] - min(spaces) + pad))
            for iratios in ratios:
                iratios[idx] = space
            # Return extra space as sum of first 2 and last 2 columns
            # of ratios from the main 5x5 gridspec
            if side in 'lt':
                return sum(ratios[0][:2])
            else:
                return sum(ratios[0][-2:])

    @_counter
    def _smart_tight_layout(self, renderer=None,
        borders=True, subplots=True, panels=True):
        """Applies special tight layout that permits flexible figure
        dimensions and preserves panel widths and subplot aspect ratios.
        The `renderer` should be a `~matplotlib.backend_bases.RendererBase`
        instance. If ``None``, renderer is inferred from
        `~matplotlib.figure.Figure.canvas.get_renderer`."""
        # Initial stuff
        axs = self._main_axes
        gridspec = self._main_gridspec
        subplots_kw = self._subplots_kw
        if subplots_kw is None or gridspec is None or not any(axs):
            return
        borders = self._tight_borders and borders
        subplots = self._tight_subplots and subplots
        panels = self._tight_panels and panels and (
            any(panel for ax in axs for panel in (ax.leftpanel,
                ax.bottompanel, ax.rightpanel, ax.toppanel))
            or any(panel for panel in (self.leftpanel,
                   self.bottompanel, self.rightpanel, self.toppanel))
            )
        if not (borders or subplots or panels):
            return

        # Tight box *around* figure
        if borders:
            # Get old un-updated bounding box
            pad = self._border_pad
            obbox = self.bbox_inches # original bbox
            oxmax, oymax = obbox.xmax, obbox.ymax
            # Get new bounding box
            bbox = self.get_tightbbox(renderer)
            xmin, xmax, ymin, ymax = bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax
            # Calculate new edges
            # TODO: Account for bounding box NaNs?
            loff, boff = xmin, ymin # left bottom margins
            roff, toff = oxmax - xmax, oymax - ymax # top right margin *deltas*
            for key,off in zip(('left','right','bottom','top'), (loff,roff,boff,toff)):
                margin = subplots_kw[key] - off + pad
                if margin < 0:
                    warnings.warn(f'Got negative {key} margin in smart tight layout.')
                elif np.isnan(margin):
                    warnings.warn('Bounding box has NaNs, cannot get outer tight layout.')
                else:
                    subplots_kw[key] = margin
            self._tight_borders_init = False

        # Prevent overlapping axis tick labels and whatnot *within* figure
        if subplots:
            # Get *original* "wspace" and "hspace" so we can adjust
            # TODO: Why not always draw figure panels, just make the gridspec
            # zero width if no panels present? Same as with axes panels.
            pad = self._subplot_pad
            olspace, orspace = subplots_kw['lspace'], subplots_kw['rspace']
            otspace, obspace = subplots_kw['tspace'], subplots_kw['bspace']
            owspace, ohspace = subplots_kw['wspace'], subplots_kw['hspace'] # originals
            owspace = [olspace, *owspace, orspace]
            ohspace = [otspace, *ohspace, obspace]
            # Get length ncols-1, nrows-1 lists, containing lists of length-2
            # lists. Each length-2 list contains groups of axes on the
            # top/bottom, left/right side of an interface that are abutt
            # against each other. For example, in 2x2 matrix, xspace will
            # contain 1 sublist with 2 length-2 lists: 1 for axes in the top
            # row, another for axes in the bottom row.
            xgroups = self._get_space_axes('x')
            ygroups = self._get_space_axes('y')
            wspace, hspace = [], []
            for ospace,groups in zip(owspace,xgroups):
                seps = [] # will use *maximum* necessary separation
                for group in groups: # groups with touching edges
                    left = max(_xtight_full(ax)[1] for ax in group[0]) # right bounds for subplot on left side of column
                    right = min(_xtight_full(ax)[0] for ax in group[1]) # left bounds for subplot on right side of column
                    seps.append((right - left)/self.dpi)
                space = max(0, ospace - min(seps) + pad) if seps else ospace
                wspace.append(space)
            for ospace, groups in zip(ohspace, ygroups):
                seps = [] # will use *maximum* necessary separation
                for group in groups:
                    top = max(_ytight_full(ax)[1] for ax in group[0])
                    bottom = min(_ytight_full(ax)[0] for ax in group[1])
                    seps.append((bottom - top)/self.dpi)
                space = max(0, ospace - min(seps) + pad) if seps else ospace
                hspace.append(space)
            # If had figure panels, need to pull out args
            # If wflush or hflush enabled, overwrite the resulting wspace and hspace
            lspace, *wspace, rspace = wspace # separate out panel spaces again
            tspace, *hspace, bspace = hspace
            if self._wflush:
                wspace = [0]*len(wspace)
            if self._hflush:
                hspace = [0]*len(hspace)
            subplots_kw.update({
                'wspace':wspace, 'hspace':hspace,
                'lspace':lspace, 'rspace':rspace, 'bspace':bspace, 'tspace':tspace,
                })

        # The same, but for spaces between figure and axes panels
        if panels:
            # Figure panels
            # Requires us to turn 2D panel axes_grids into lists of stacked panels
            for side in 'lrbt':
                paxs = getattr(self, side + 'panel')
                if not paxs:
                    continue
                if side in 'bt':
                    paxs = [paxs[:,i] for i in range(paxs.shape[1])]
                else:
                    paxs = [paxs[i,:] for i in range(paxs.shape[0])]
                sep = self._panel_tight_layout(side, paxs, renderer, figure=True)
                subplots_kw[side + 'sep'] = sep
            # Adjust space in axes panel wratios/hratios lists, return space
            # that must be alotted for panels in aspect ratio conservation calc
            wpanels, hpanels = [], []
            baxs = self._get_side_axes('b') # lists of baxs for each row
            taxs = self._get_side_axes('t')
            laxs = self._get_side_axes('l')
            raxs = self._get_side_axes('r')
            for ibaxs,itaxs in zip(baxs,taxs):
                ihpanels = 0
                ihpanels += self._panel_tight_layout('b', ibaxs, renderer)
                ihpanels += self._panel_tight_layout('t', itaxs, renderer)
                hpanels.append(ihpanels)
            for ilaxs,iraxs in zip(laxs,raxs):
                iwpanels = 0
                iwpanels += self._panel_tight_layout('l', ilaxs, renderer)
                iwpanels += self._panel_tight_layout('r', iraxs, renderer)
                wpanels.append(iwpanels)
            # Room for panels, ignoring the figure panel gridspec positions
            _, *wpanels, _ = wpanels
            _, *hpanels, _ = hpanels
            subplots_kw.update({'wpanels':wpanels, 'hpanels':hpanels})

        # Update gridspec(s)
        # self.subplots_adjust(**{key:value for key,value in gridspec_kw.items() if key in ('left','right','top','bottom')})
        figsize, gridspec_kw, _ = _subplots_geometry(**subplots_kw)
        gridspec.update(**gridspec_kw)
        self.set_size_inches(figsize, forward=True)
        self._update_ratios(**gridspec_kw) # update axes width/height in wratios/hratios
        # Axes width and height used for unit scaling in various places
        width, height = figsize
        for ax in axs:
            width_new = abs(ax._position.width)*width
            height_new = abs(ax._position.height)*height
            ax.width, ax.height = width_new, height_new

    def _sync_panels(self, ax, side, *, space, width, sep):
        """Syncs panel settings across rows and columns of the figure.
        Requires input width, separation, and space arrays."""
        # Axis sharing setup
        # TODO: Do this at end only? Nah, need future user actions to sync
        pax = getattr(ax, side + 'panel')
        if not pax:
            return # no panels so nothing to be done
        if pax[0]._share:
            ax._share_panels_setup()
        group = self._get_shared_axes('y' if side in 'lr' else 'x', ref=ax)
        parent = group[0]
        for child in group[1:]:
            if side in 'lr':
                child._sharey_setup(parent, ax._sharey_level)
            else:
                child._sharex_setup(parent, ax._sharex_level)

        # Modify gridspec in all panels on the group
        nrows, ncols = pax[0]._stack_gridspec.get_geometry()
        axs = self._get_side_axes(side, ref=ax)
        for iax in axs:
            # Fix main gridspec
            gridspec = iax._main_gridspec
            if side in 'lr':
                iratios = gridspec.get_width_ratios()
            else:
                iratios = gridspec.get_height_ratios()
            idx = (1 if side in 'lt' else -2)
            iratios[idx] = space # space between subplot and panels
            idx = (0 if side in 'lt' else -1)
            iratios[idx] = sum(sep) + sum(width)
            # Fix axes ratio, so it is in physical units
            # iratios[2] -= (sum(sep) + sum(width) + space)
            # Fix stack gridspec ratios
            pax = getattr(iax, side + 'panel')
            if not pax: # sync seps
                continue
            igridspec = pax[0]._stack_gridspec
            inrows, incols = igridspec.get_geometry()
            if nrows != inrows or ncols != incols:
                raise ValueError(f'Conflicting number of stacked axes or panels in this row or column.')
            if side in 'lr':
                ipratios = igridspec.get_width_ratios()
            else:
                ipratios = igridspec.get_height_ratios()
            ipratios[::2] = width
            ipratios[1::2] = sep

        # Modify wpanels and hpanels, and sync across rows and cols
        # TODO: Fix aspect ratios in light of wpanels and hpanels at
        # draw time with the aspect_fix function?
        idx = (0 if side in 'lt' else -1)
        if side in 'lr':
            key = 'wpanels'
            loc = _xrange(ax)[idx]-1 # -1 excludes fig panels
        else:
            key = 'hpanels'
            loc = _yrange(ax)[idx]-1
        # Apply changes
        # TODO: Is there ever situation where += will mess things up?
        gridspec = self._main_gridspec
        subplots_kw = self._subplots_kw
        subplots_kw[key][loc] += sum(sep) + sum(width) + space
        figsize, gridspec_kw, _ = _subplots_geometry(**subplots_kw)
        gridspec.update(**gridspec_kw)
        self.set_size_inches(figsize, forward=True)
        self._update_ratios(**gridspec_kw)

    def _update_aspect(self):
        """Adjust average aspect ratio used for gridspec calculations. This
        fixes grids with identically fixed aspect ratios, e.g. identically
        zoomed-in cartopy projections and imshow images."""
        if not self._main_axes:
            return
        ax = self._main_axes[self._ref_num-1]
        aspect = None
        subplots_kw = self._subplots_kw
        if isinstance(ax, axes.CartopyProjectionAxes):
            # Basemap projections are static but cartopy can change
            bbox = ax.background_patch._path.get_extents()
            aspect = abs(bbox.width)/abs(bbox.height)
        elif isinstance(ax, axes.CartesianAxes):
            # Image plots, etc. Mostly copied from apply_aspect().
            mode = ax.get_aspect()
            xscale, yscale = ax.get_xscale(), ax.get_yscale()
            if mode == 'equal':
                if xscale == 'log' and yscale == 'log':
                    aspect = 1.0/ax.get_data_ratio_log()
                elif xscale == 'linear' and yscale == 'linear':
                    aspect = 1.0/ax.get_data_ratio()
                else:
                    pass # matplotlib issues warning, forces aspect == 'auto'
        # Apply aspect
        if aspect is not None and aspect != subplots_kw['aspect']:
            gridspec = self._main_gridspec
            subplots_kw['aspect'] = aspect
            figsize, gridspec_kw, _ = _subplots_geometry(**subplots_kw)
            gridspec.update(**gridspec_kw)
            self.set_size_inches(figsize)
            self._update_ratios(**gridspec_kw)

    def _update_axislabels(self, axis=None, span=False, **kwargs):
        """Get axis label for axes with axis sharing or spanning enabled.
        When `span` is False, we add labels to all axes. When `span` is
        True, we filter to one axes."""
        # Get spanning axes
        # TODO: Use fig.align_ylabels to make sure global y label is *always*
        # properly offset? Panels may make this too complicated.
        # TODO: Aligned axis labels for arbitrary sides, e.g. right labels
        # and top labels! Is this already a thing!?
        axes = [axis]
        record = self._axes_spanning
        if axis is None:
            axes = record
        # Loop through axes
        for axis in axes:
            base = axis.axes
            name = axis.axis_name
            side = ('b' if name == 'x' else 'l')
            along = (_xrange if name=='x' else _yrange)
            if name not in 'xy': # i.e. theta or radius
                return
            for _ in range(2): # try 2 levels down, should be sufficient
                base = getattr(base, '_share' + name, None) or base
            # Get the x/y axes, make sure to select shared panels if possible
            spanon = getattr(base, '_span'  + name)
            if not spanon:
                axes = [getattr(base, name + 'axis')]
            else:
                axs = self._get_side_axes(side, ref=base)
                axs = [getattr(ax, '_share' + name) or ax for ax in axs]
                axes = [getattr(ax, name + 'axis') for ax in axs]
            # Apply settings and store list of spanning axes, and change
            # vertical align from baseline to bottom so less cramped
            for axis in axes:
                if spanon and axis not in record:
                    record.append(axis)
                axis.label.set_visible(True)
                axis.label.update(kwargs)
            # If requested, turn on spanning
            # TODO: Use tightbbox for this?
            if not span or not spanon:
                continue
            # Use halfway-point axis for spanning label, because for odd
            # number of columns or rows, gives us perfect offset
            ranges = np.array([along(ax) for ax in axs])
            saxis = axes[(np.argmin(ranges[:,0]) + np.argmax(ranges[:,0]))//2]
            for axis in axes:
                if saxis is not axis:
                    axis.label.set_visible(False)
            # Reposition to "span" other axes
            idx = slice(ranges.min(), ranges.max() + 1)
            if name == 'x': # span columns
                subspec = self._main_gridspec[0,idx]
            else: # spans rows
                subspec = self._main_gridspec[idx,0]
            bbox = subspec.get_position(self) # in figure-relative coordinates
            x0, y0, width, height = bbox.bounds
            if name == 'x':
                transform = mtransforms.blended_transform_factory(self.transFigure, mtransforms.IdentityTransform())
                position = (x0 + width/2, 1)
            else:
                transform = mtransforms.blended_transform_factory(mtransforms.IdentityTransform(), self.transFigure)
                position = (1, y0 + height/2)
            saxis.label.update({'position':position, 'transform':transform})

    def _update_ratios(self, **kwargs):
        """Updates the width and height ratios values in the 5x5 subplot
        gridspecs to account for changes to figure size -- makes sure
        ratio values match physical units in inches. Without this, if the
        figure size changes, panel widths change! Keyword args are sent
        to `~matplotlib.gridspec.GridSpec.update`."""
        # Main gridspec ratios
        gridspec = self._main_gridspec
        wratios = gridspec.get_width_ratios()
        hratios = gridspec.get_height_ratios()
        for ax in self._main_axes:
            # Get extra space for panels
            igridspec = ax._main_gridspec
            iwratios = igridspec.get_width_ratios() # gets my *custom* ratios
            ihratios = igridspec.get_height_ratios()
            wpanels = sum(r for i,r in enumerate(iwratios) if i != 2)
            hpanels = sum(r for i,r in enumerate(ihratios) if i != 2) # all entries of 5x5 gridspec
            # Get current space alotted in main gridspec
            xrange = 2*np.array(_xrange(ax)) # get range for *actual* grid, not just visible grid
            yrange = 2*np.array(_yrange(ax))
            fullwidth = sum(wratios[xrange[0]:xrange[1]+1])
            fullheight = sum(hratios[yrange[0]:yrange[1]+1])
            # Adjust child gridspec to reflect changes in parent
            iwratios[2] = fullwidth - wpanels
            ihratios[2] = fullheight - hpanels
            igridspec.set_width_ratios(iwratios)
            igridspec.set_height_ratios(ihratios)
        # Update main gridspec to reflect changes to child gridspecs.
        # Required because update propagates *down*, not *up*
        # NOTE: gridspec.update *always* needs the original keyword args, or
        # it seems to impose default properties!
        gridspec.update(**kwargs)

    def _update_suplabels(self, ax, labels, rows=True, **kwargs):
        """Assigns row labels and column labels, updates label settings."""
        attr = ('rowlabel' if rows else 'collabel')
        _along = (_yrange if rows else _xrange) # along labels
        _across = (_xrange if rows else _yrange) # across from labels
        if not self._main_axes:
            axs = [ax]
        else:
            axs = [ax for ax in self._main_axes if _across(ax)[0] == 1] # axes on the edge
            axs = [ax for _,ax in sorted(zip([_along(ax)[0] for ax in axs], axs))] # order by yrange
        if labels is None or isinstance(labels, str): # common during testing
            labels = [labels]*len(axs)
        if len(labels) != len(axs):
            raise ValueError(f'Got {len(labels)} {attr}s, but there are {len(axs)} {attr.split("label")[0]}s.')
        for ax,label in zip(axs,labels):
            obj = getattr(ax, attr)
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

    def add_subplot(self, *args, **kwargs):
        """Establishes `~proplot.axes.CartesianAxes` as the default projection.
        Issues warning for new users that try to access the
        `~matplotlib.figure.Figure.add_subplot` and
        `~matplotlib.figure.Figure.colorbar` functions."""
        kwargs.setdefault('projection', 'cartesian')
        if self._locked:
            warnings.warn('Using "add_subplot" or "colorbar" with ProPlot figures may result in unexpected behavior. '
                'Please use the subplots() command to create your figure, subplots, and panels all at once.')
        ax = super().add_subplot(*args, **kwargs)
        return ax

    def draw(self, renderer):
        """Applies various post-processing steps, including tight layout
        adjustment and aspect ratio preservation."""
        if self._post_init:
            self._update_aspect()
            self._update_layout(renderer) # want user to have ability to call it manually
            self._update_axislabels(span=True)
        super().draw(renderer)

    def savefig(self, filename, **kwargs):
        """
        Applies various post-processing steps, including tight layout
        adjustment and aspect ratio preservation. Then calls the parent
        `~matplotlib.figure.Figure.savefig` method.

        Parameters
        ----------
        filename : str
            The file path. User directories are automatically
            expanded, e.g. ``fig.save('~/plots/plot.png')``.
        **kwargs
            Passed to `~matplotlib.figure.Figure.savefig`.
        """
        filename = os.path.expanduser(filename)
        if self._post_init:
            self._update_aspect()
            self._update_layout() # necessary! get weird layout without this
            self._update_axislabels(span=True)
        super().savefig(filename, **kwargs)

    save = savefig
    """Alias for `~Figure.savefig`... because calling ``fig.savefig``
    is sort of redundant."""

    # Define panels as immutable properties, can only be set internally
    # This also documents the properties in sphinx
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


#-----------------------------------------------------------------------------#
# Primary plotting function, used to create figure/axes
#-----------------------------------------------------------------------------#
def _panels_kwargs(
    panels, colorbars, legends,
    panels_kw, colorbars_kw=None, legends_kw=None,
    figure=False,
    ):
    """Returns dictionary with keyword `[lrbt]space`, `[lrbt]width`,
    `[lrbt]sep`, `[lrbt]flush`, `[lrbt]share`, and for figure panels,
    `[lrbt]array`, with defaults filled in."""
    # Get which panels
    kwargs = {}
    panels = _side_translate_inverse.get(panels, panels)
    legends = _side_translate_inverse.get(legends, legends)
    colorbars = _side_translate_inverse.get(colorbars, colorbars)
    onpanels = panels + colorbars + legends
    allpanels = 'lrbt'
    if not {*onpanels} <= {*allpanels}:
        raise ValueError(f'Invalid panel spec "{onpanels}" for {"figure" if figure else "axes"} panels. Valid characters are: {", ".join(allpanels)}.')
    if len({*onpanels}) != len(onpanels):
        raise ValueError('You requested the same panel side more than once, e.g. by passing panels="b" and colorbars="b".')
    if colorbars_kw is None: # they are empty dict when this function is used to build axes panels
        colorbars_kw = panels_kw
    if legends_kw is None:
        legends_kw = panels_kw

    # Fill non-panels with empty args, to populate dict with args required
    # by _subplots_geometry
    names = ('space', 'width', 'sep', 'flush', 'share')
    if figure:
        names = (*names, 'array')
    for side in {*allpanels} - {*onpanels}: # defaults for *absent* panels!
        for name in names:
            kwargs[side + name] = None
    # Add unknown arguments to output, will raise errors down the line
    names = (*names, 'stack') # just used to set the 'sep' and 'width' arrays
    regex = re.compile(f'^[tlrb]?({"|".join(names)})')
    for kw in (panels_kw, colorbars_kw, legends_kw):
        for key,value in kw.items():
            if not regex.match(key):
                kwargs[key] = value

    # Helper func that returns property from the appropriate dictionary, or the
    # default from rc.subplots.
    def _prop(side, key, defaults):
        if not isinstance(defaults, tuple):
            defaults = 3*(defaults,)
        for which,kw,default in zip((panels, colorbars, legends), (panels_kw, colorbars_kw, legends_kw), defaults):
            if side not in which:
                continue
            if isinstance(default, str):
                default = units(rc['subplots.' + default])
            vglobal = kw.get(key, None)
            vside = kw.get(side + key, None)
            if vglobal is not None and vside is not None:
                warnings.warn(f'Got conflicting values {key}={vglobal} and {side}{key}={vside}. Using "{key}".')
            return _default(vglobal, vside, default)
    # Apply default props for all panels that are turned on
    for side in onpanels:
        # Axis sharing
        share = _prop(side, 'share', (True,False,False))
        kwargs[side + 'share'] = share
        # Widths
        width = _prop(side, 'width', ('panelwidth', 'cbarwidth', 'legwidth'))
        width = np.atleast_1d(units(width))
        stack = _prop(side, 'stack', len(width) if np.iterable(width) else stack)
        if stack < 1:
            raise ValueError(f'"{side+stack}" argument must be integer >=1.')
        if len(width) == 1:
            width = np.repeat(width, (stack,))
        if len(width) != stack:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(width)} widths.')
        kwargs[side + 'width'] = width
        # Figure or axes panel-specific props
        if figure:
            kwargs[side + 'array'] = _prop(side, 'array', None)
            space = ('xlabspace' if side == 'b' else 'ylabspace' if side == 'l' else 'nolabspace')
        else:
            space = ('panelspace' if share else 'xlabspace' if side == 'b' else 'ylabspace' if side == 'l' else 'panelspace')
        kwargs[side + 'space'] = units(_prop(side, 'space', space))
        # Stacked panels
        flush = _prop(side, 'flush', False)
        if stack == 1 or flush:
            sep = 0
        else:
            default = 'nolabspace' if share else 'ylabspace' if side in 'lr' else 'xlabspace'
            sep = _prop(side, 'sep', (default,default,default))
        sep = np.atleast_1d(units(sep))
        if len(sep) == 1:
            sep = np.repeat(sep, (stack-1,))
        if len(sep) != stack-1:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(sep)} separations.')
        kwargs[side + 'sep'] = sep
        kwargs[side + 'flush'] = flush
    if not figure:
        kwargs['sides'] = onpanels
    return kwargs

def _subplots_geometry(**kwargs):
    """Saves arguments passed to `subplots`, calculates gridspec settings and
    figure size necessary for requested geometry, returns keyword args
    necessary to reconstruct and modify this configuration, as is done during
    tight layout scaling."""
    # Necessary arguments to reconstruct this grid, with defaults filled in
    # NOTE: Why put everything into kwargs dictionary? Because need to store
    # the dictionary for tight layout scaling.
    subplots_kw = {}
    def _pop(key): # get from array and apply defaults
        if key not in kwargs:
            raise ValueError(f'Required argument "{key}" not passed to _subplots_geometry.')
        value = kwargs.pop(key)
        subplots_kw[key] = value # record for future calls to _subplots_geometry
        return value
    # Basic dimensions and geometry
    nrows, ncols       = _pop('nrows'), _pop('ncols')
    aspect, xref, yref = _pop('aspect'), _pop('xref'), _pop('yref')
    width, height      = _pop('width'), _pop('height')
    axwidth, axheight  = _pop('axwidth'), _pop('axheight')
    # Space between subplots and space inside each row/column allocated
    # for axes panels.
    wpanels, hpanels = _pop('wpanels'), _pop('hpanels')
    hspace, wspace   = _pop('hspace'), _pop('wspace')
    hratios, wratios = _pop('hratios'), _pop('wratios')
    left, bottom     = _pop('left'), _pop('bottom')
    right, top       = _pop('right'), _pop('top')
    # Various panel settings
    # Some are only needed when setting up panel axes, so we ignore them
    _, bwidth, bspace = _pop('barray'), _pop('bwidth'), _pop('bspace')
    _, lwidth, lspace = _pop('larray'), _pop('lwidth'), _pop('lspace')
    _, rwidth, rspace = _pop('rarray'), _pop('rwidth'), _pop('rspace')
    _, twidth, tspace = _pop('tarray'), _pop('twidth'), _pop('tspace')
    rsep, _, _ = _pop('rsep'), _pop('rshare'), _pop('rflush')
    bsep, _, _ = _pop('bsep'), _pop('bshare'), _pop('bflush')
    lsep, _, _ = _pop('lsep'), _pop('lshare'), _pop('lflush')
    tsep, _, _ = _pop('tsep'), _pop('tshare'), _pop('tflush')
    # Defaults for below calculations
    lsep, rsep = _default(lsep, []), _default(rsep, [])
    tsep, bsep = _default(tsep, []), _default(bsep, [])
    lspace, rspace = _default(lspace, 0), _default(rspace, 0)
    tspace, bspace = _default(tspace, 0), _default(bspace, 0)
    lwidth, rwidth = _default(lwidth, [0]), _default(rwidth, [0])
    twidth, bwidth = _default(twidth, [0]), _default(bwidth, [0])
    # Check no extra args passed
    if kwargs:
        raise ValueError(f'Unknown keyword args passed to _subplots_geometry: {kwargs}.')

    # Reference properties
    # NOTE: Axes range is for gridspec including panels, but wratios and
    # wspace arrays do not include panels. Below accounts for this. Don't
    # want to add e.g. lwidth to wratios yet because don't yet know physical
    # space that wratios subplots will occupy.
    dx, dy = xref[1]-xref[0], yref[1]-yref[0]
    rwspace = sum(wspace[xref[0]:xref[1]-1])
    rhspace = sum(hspace[yref[0]:yref[1]-1])
    rwratio = ncols*(sum(wratios[xref[0]:xref[1]])/dx)/sum(wratios)
    rhratio = nrows*(sum(hratios[yref[0]:yref[1]])/dy)/sum(hratios)
    bpanel = sum(bwidth) + sum(bsep) + bspace
    rpanel = sum(rwidth) + sum(rsep) + rspace
    lpanel = sum(lwidth) + sum(lsep) + lspace
    tpanel = sum(twidth) + sum(tsep) + tspace
    if np.iterable(aspect):
        aspect = aspect[0]/aspect[1]

    # Determine figure dims from axes width and height, or vice versa
    # NOTE: Account for 'whitespace' that is spanned by our reference axes.
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    if width is None and height is None: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = units(rc['subplots.axwidth'])
        if axheight is not None:
            auto_width = True
            axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows
            height = axheight_all + top + bottom + sum(hspace) + sum(hpanels) + bpanel + tpanel
        if axwidth is not None:
            auto_height = True
            axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
            width = axwidth_all + left + right + sum(wspace) + sum(wpanels) + rpanel + lpanel
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
    else:
        if height is not None:
            axheight_all = height - top - bottom - sum(hspace) - sum(hpanels) - bpanel - tpanel
            axheight = (axheight_all*dy/nrows) + rhspace # reverse engineered from above
        if width is not None:
            axwidth_all = width - left - right - sum(wspace) - sum(wpanels) - rpanel - lpanel
            axwidth = (axwidth_all*dx/ncols) + rwspace

    # Automatically scale fig dimensions
    # For e.g. common use case [[1,1,2,2],[0,3,3,0]], make sure we still scale
    # the reference axes like a square even though it occupes two columns of gridspec!
    if auto_width:
        axwidth = axheight*aspect
        axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
        width = axwidth_all + left + right + sum(wspace) + sum(wpanels) + rpanel + lpanel
    elif auto_height:
        axheight = axwidth/aspect
        axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows
        height = axheight_all + top + bottom + sum(hspace) + sum(hpanels) + bpanel
    if axwidth_all < 0:
        raise ValueError(f"Not enough room for axes (would have width {axwidth_all}). Try using tight=False, increasing figure width, or decreasing 'left', 'right', or 'wspace' spaces.")
    if axheight_all < 0:
        raise ValueError(f"Not enough room for axes (would have height {axheight_all}). Try using tight=False, increasing figure height, or decreasing 'top', 'bottom', or 'hspace' spaces.")

    # Keyword args for gridspec class
    # Make 'ratios' and 'spaces' in physical units
    wspace = [lspace, *wspace, rspace] # may be zero
    hspace = [tspace, *hspace, bspace]
    lpanel, rpanel = sum(lwidth) + sum(lsep), sum(rwidth) + sum(rsep)
    tpanel, bpanel = sum(twidth) + sum(tsep), sum(bwidth) + sum(bsep)
    wratios = [lpanel, *(axwidth_all*wratios/sum(wratios) + wpanels), rpanel]
    hratios = [tpanel, *(axheight_all*hratios/sum(hratios) + hpanels), bpanel]
    left, right = left/width, 1 - right/width
    top, bottom = 1 - top/height, bottom/height
    gridspec_kw = {
        'nrows': nrows+2, 'ncols': ncols+2, # plus 2 is room for panels
        'left': left, 'bottom': bottom, 'right': right, 'top': top,
        'wspace': wspace, 'hspace': hspace,
        'width_ratios': wratios, 'height_ratios' : hratios,
        }
    return (width, height), gridspec_kw, subplots_kw

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

def _args_string(name, args):
    """String for ignored arguments warning."""
    return f'You passed {name}=True with ' \
        + ', '.join(f'{key}={value}' for key,value in args.items()) + '. ' \
        + ', '.join(f'"{key}"' for key in args.keys()) + ' will be ignored.'

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
        raise ValueError(f'Unknown journal figure size specifier "{journal}". ' +
                          'Current options are: ' + ', '.join(table.keys()))
    # Return width, and optionally also the height
    width, height = None, None
    try:
        width, height = value
    except TypeError:
        width = value
    return width, height

def subplots(array=None, ncols=1, nrows=1,
    ref=1, # reference axes for fixing aspect ratio
    order='C', # allow calling with subplots(array)
    aspect=1, figsize=None,
    width=None,  height=None, axwidth=None, axheight=None, journal=None,
    wwidth=None, hwidth=None,
    axwidths=None, axheights=None,
    hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
    width_ratios=None, height_ratios=None,
    flush=None, wflush=None, hflush=None,
    left=None, bottom=None, right=None, top=None, # spaces around edge of main plotting area, in inches
    tight=None, tightborders=None, tightsubplots=None, tightpanels=None,
    borderpad=None, panelpad=None, subplotpad=None,
    span=None, spanx=None, spany=None, # custom setting, optionally share axis labels for axes with same xmin/ymin extents
    share=None, sharex=None, sharey=None, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
    panel=None, panels=None, legend=None, legends=None, colorbar=None, colorbars=None,
    axpanel=None, axlegend=None, axcolorbar=None,
    axpanel_kw=None, axcolorbar_kw=None, axlegend_kw=None,
    axpanels=None, axlegends=None, axcolorbars=None,
    axpanels_kw=None, axcolorbars_kw=None, axlegends_kw=None,
    basemap=False, proj=None, projection=None, proj_kw=None, projection_kw=None,
    autoformat=True, # arguments for figure instantiation
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
        the order axes appear in the `axs` list, and the order of a-b-c
        labelling when using `~proplot.axes.BaseAxes.format` with ``abc=True``.
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
    hratios, wratios, axheights, axwidths : optional
        Aliases for `height_ratios`, `width_ratios`.
    height_ratios, width_ratios : float or list thereof, optional
        Passed to `FlexibleGridSpecBase`. The height
        and width ratios for the subplot grid. Length of `height_ratios`
        must match the number of rows, and length of `width_ratios` must
        match the number of columns.
    hspace, wspace : float or str or list thereof, optional
        If passed, turns off `tightsubplots`.
        These are passed to `FlexibleGridSpecBase`, denote the
        spacing between each column and row of the grid. If float
        or string, expanded into lists of length ``ncols-1`` (for `wspace`) or
        length ``nrows-1`` (for `hspace`). For each element of the list, if float,
        units are inches. If string, units are interpreted by `~proplot.utils.units`.
    top, bottom, left, right : float or str, optional
        If passed, turns off `tightborders`. These are passed to
        `FlexibleGridSpecBase`, denote the width of padding between the subplots
        and the figure edge. If float, units are inches. If string, units are
        interpreted by `~proplot.utils.units`.

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

    proj, projection : str or dict-like, optional
        The map projection name. If string, applies to all subplots. If
        dictionary, values apply to specific subplots, as with `axpanels`.
        If an axes projection is not specified in the dictionary, that axes
        will be Cartesian.

        For example, with ``ncols=4`` and ``proj={1:'mercator', (2,3):'hammer'})``,
        the leftmost subplot is a Mercator projection, the middle 2 are
        Hammer projections, and the rightmost is a normal Cartesian axes.
    proj_kw, projection_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class on instantiation.
        If dictionary of properties, applies globally. If *dictionary of
        dictionaries* of properties, applies to specific subplots, as with `axpanels`.

        For example, with ``ncols=2`` and ``proj_kw={1:{'lon_0':0}, 2:{'lon_0':180}}``,
        the projection in the left subplot is centered on the prime meridian,
        and the projection in the right subplot is centered on the
        international dateline.
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.
        If boolean, applies to all subplots. If dictionary, values apply to
        specific subplots, as with `axpanels`.

    panel : str, optional
        Specify which sides of the figure should have a "panel".
        This is great for creating global colorbars and legends.
        String should contain any of the characters ``'l'`` (left panel),
        ``'r'`` (right panel), ``'b'`` (bottom panel), or ``'t'`` (top panel).
        For example, ``'br'`` will draw a right and bottom panel.

        Panel axes are stored as ``leftpanel``, ``rightpanel``,
        ``bottompanel``, and ``toppanel`` attributes on the figure object.
        They can also be accessed by the attribute aliases ``lpanel``,
        ``rpanel``, ``bpanel``, and ``tpanel``.
    panels : str, optional
        As with `panel`, but the default behavior is to assign a panel
        to *every* row or column of subplots. Individual panels can then
        be accessed with e.g. ``fig.leftpanel[0]``, ``fig.leftpanel[1]``.
    colorbar, legend, colorbars, legends : optional
        Identical to `panel` and `panels`, except the *default* panel width is
        more appropriate for being "filled" with a colorbar or legend with
        the `~proplot.axes.PanelAxes.colorbar` and
        the `~proplot.axes.PanelAxes.legend` methods.
    larray, rarray, barray, tarray : list of int, optional
        Defines how figure panels span rows and columns of subplots.
        Interpreted like `array` -- the integers specify panels that span
        *arbitrary, contiguous* columns or rows of subplots.

        For example, ``plot.suplots(ncols=3, panels='b', barray=[1, 2, 2])``
        draws a panel on the bottom of the first column and spanning the bottom
        of the right 2 columns, and ``barray=[0, 2, 2]`` only draws a panel
        underneath the right 2 columns -- as with `array`, the ``0`` indicates
        an empty space.
    lwidth, rwidth, bwidth, twidth : optional
        See `~proplot.axes.BaseAxes.panel_axes`, usage is identical.
    lstack, rstack, bstack, tstack : optional
        See `~proplot.axes.BaseAxes.panel_axes`, usage is identical.
    lsep, rsep, bsep, tsep : optional
        See `~proplot.axes.BaseAxes.panel_axes`, usage is identical.
    lshare, rshare, bshare, tshare : optional
        See `~proplot.axes.BaseAxes.panel_axes`, usage is identical.
    lspace, rspace, bspace, tspace : float, optional
        If passed, turns off `tightsubplots`. As in
        `~proplot.axes.BaseAxes.panel_axes`, but controls space between the
        main subplot grid and the figure panels.
    lflush, rflush, bflush, tflush : optional
        As in `~proplot.axes.BaseAxes.panel_axes`, but only controls whether
        *stacked* panels are flush against each other -- i.e. does not make
        figure panels flush against the main subplots.

    axpanel, axpanels : str or dict-like, optional
        Bulk adds axes panels with `~proplot.axes.BaseAxes.panel_axes`. You
        can also build on-the-fly panels in myriad ways (see
        :ref:`On-the-fly panels`). Both `axpanel` and `axpanels`
        are acceptable. The argument is interpreted as follows:

        * If string, panels are drawn on the same side for all subplots.
          String should contain any of the characters ``'l'`` (left panel),
          ``'r'`` (right panel), ``'t'`` (top panel), or ``'b'`` (bottom panel).
          For example, ``'rt'`` will draw a right and top panel.
        * If dict-like, panels can be drawn on different sides for
          different subplots. For example, for a 4-subplot figure,
          ``axpanels={1:'r', (2,3):'l'}`` indicates that we want to
          draw a panel on the right side of subplot number 1, on the left
          side of subplots 2 and 3, and **no panel** on subplot 4.

        Panel axes are stored as ``leftpanel``, ``rightpanel``,
        ``bottompanel``, and ``toppanel`` attributes on axes objects.
        They can also be accessed by the attribute aliases ``lpanel``,
        ``rpanel``, ``bpanel``, and ``tpanel``.
    axcolorbar, axcolorbars : optional
        Identical to `axpanel` and `axpanels`, except the *default* panel
        width is more appropriate for a colorbar. The panel can then be
        "filled" with a colorbar with e.g. ``ax.rpanel.colorbar()``.
    axlegend, axlegends : optional
        Identical to `axpanel` and `axpanels`, except the *default* panel
        width is more appropriate for a legend. The panel can then be
        "filled" with a legend with e.g. ``ax.rpanel.legend()``.
    axpanel_kw, axpanels_kw : dict-like, optional
        Keyword args passed to `~proplot.axes.BaseAxes.panel_axes` for panels
        listed in the `axpanel` and `axpanels` keyword args.
        If dictionary of properties, applies globally. If *dictionary of
        dictionary* of properties, applies to specific subplots, as with `axpanels`.

        For example, consider a 2-subplot figure with ``axpanels='l'``.
        With ``{'lwidth':1}``, both left panels will be 1 inch wide.
        With ``{1:{'lwidth':1}, 2:{'lwidth':0.5}}``, the left subplot
        panel will be 1 inch wide and the right subplot panel will be
        0.5 inches wide.
    axcolorbar_kw, axcolorbars_kw : optional
        As with `axpanel_kw`, but for panels listed in the `axcolorbar`
        and `axcolorbars` keyword args.
    axlegend_kw, axlegends_kw : optional
        As with `axpanel_kw`, but for panels listed in the `axlegend`
        and `axlegends` keyword args.

    Other parameters
    ----------------
    tight, tightborders, tightsubplots, tightpanels, borderpad, subplotpad, panelpad, flush, wflush, hflush, autoformat
        Passed to `Figure`. Defaults are as follows.

        * `tightborders` defaults to ``False`` if user provided the `top`,
          `bottom`, `left`, or `right` keyword args. ``True`` otherwise.
        * `tightsubplots` defaults to ``False`` if user  provided the `wspace`
          or `hspace` gridspec keyword args, or the `lspace`, `rspace`, `bspace`,
          `lsep`, `rsep`, or `bsep` figure panel keyword args. ``True`` otherwise.
        * `tightpanels` defaults to ``False`` if user provided the `lspace`,
          `rspace`, `tspace`, `bspace`, `lsep`, `rsep`, `tsep`, or `bsep` axes
          panel keyword args with e.g. the `axpanels_kw` dictionary. ``True``
          otherwise.

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
    # TODO: Generalize axes sharing for right y-axes and top x-axes. Enable a secondary
    # axes sharing mode where we *disable ticklabels and labels*, but *do not
    # use the builtin sharex/sharey API*, suitable for complex map projections.
    # For spanning axes labels, right now only detect **x labels on bottom**
    # and **ylabels on top**. Generalize for all subplot edges.
    rc._getitem_mode = 0 # ensure still zero; might be non-zero if had error in 'with context' block
    #-------------------------------------------------------------------------#
    # Array setup
    #-------------------------------------------------------------------------#
    # Build array
    if order not in ('C','F'): # better error message
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
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
    # Axes panels
    #-------------------------------------------------------------------------#
    # Flexible args
    axpanels    = _default(axpanel, axpanels, '')
    axcolorbars = _default(axcolorbar, axcolorbars, '')
    axlegends   = _default(axlegend, axlegends, '')
    axpanels_kw    = _default(axpanel_kw, axpanels_kw, {})
    axcolorbars_kw = _default(axcolorbar_kw, axcolorbars_kw, {})
    axlegends_kw   = _default(axlegend_kw, axlegends_kw, {})
    for iname,ipanels,ikw in (
        ('panels', axpanels, axpanels_kw),
        ('colorbars', axcolorbars, axcolorbars_kw),
        ('legends', axlegends, axlegends_kw)
        ):
        if not ipanels and ikw:
            warnings.warn(f'ax{iname}={repr(ipanels)}, ignoring ax{iname}_kw: {ikw}.')
    # Create dictionaries of panel toggles and settings
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    axpanels    = _axes_dict(naxs, axpanels, kw=False, default='')
    axcolorbars = _axes_dict(naxs, axcolorbars, kw=False, default='')
    axlegends   = _axes_dict(naxs, axlegends, kw=False, default='')
    axpanels_kw    = _axes_dict(naxs, axpanels_kw, kw=True)
    axcolorbars_kw = _axes_dict(naxs, axcolorbars_kw, kw=True)
    axlegends_kw   = _axes_dict(naxs, axlegends_kw, kw=True)
    # Turn off tightborders
    args = {key:value for key,value in
        (('left',left),('right',right),('top',top),('bottom',bottom))
        if value is not None}
    if args:
        tightborders = _default(tightborders, False)
        if tightborders:
            warnings.warn(_args_string('tightborders', args))
    # Turn off tightsubplots
    args = {**{key:value for key,value in (('wspace',wspace),('hspace',hspace))
        if value is not None}, **{key:value for key,value in kwargs.items()
        if ('space' in key or 'sep' in key) and value is not None
        and not (np.iterable(value) and not len(value))}}
    if args:
        tightsubplots = _default(tightsubplots, False)
        if tightsubplots:
            warnings.warn(_args_string('tightsubplots', args))
    # Turn off tightpanels
    args = {key:value for ikw in (axpanels_kw,axcolorbars_kw,axlegends_kw)
        for jkw in ikw.values() for key,value in jkw.items()
        if ('space' in key or 'sep' in key) and value is not None
        and not (np.iterable(value) and not len(value))}
    if args:
        tightpanels = _default(tightpanels, False)
        if tightpanels:
            warnings.warn(_args_string('tightpanels', args))

    #-------------------------------------------------------------------------#
    # Figure and axes panels
    #-------------------------------------------------------------------------#
    # Get axes panel keyword args
    for num in range(1,naxs+1):
        axpanels_kw[num] = _panels_kwargs(
            axpanels[num], axcolorbars[num], axlegends[num],
            axpanels_kw[num], axcolorbars_kw[num], axlegends_kw[num],
            figure=False)
    # Get default figure panel keyword args
    panels = _default(panel, panels, '')
    legends = _default(legend, legends, '')
    colorbars = _default(colorbar, colorbars, '')
    kwargs = _panels_kwargs(panels, colorbars, legends, kwargs, figure=True)
    # Get figure panel arrays
    for ipanel,ispan in zip((panel,legend,colorbar,panels,legends,colorbars),(1,1,1,0,0,0)):
        if ipanel is not None and not isinstance(ipanel, str):
            raise ValueError(f'Figure panel input must be string containing any of the characters: l, r, b.')
        for side in (ipanel or ''):
            value = kwargs.get(side + 'array', None)
            nmax, name = ((ncols,'cols') if side in 'bt' else (nrows,'rows'))
            if value is None:
                if ispan:
                    value = [1]*nmax
                else:
                    value = [*range(1,nmax+1)]
            elif not np.iterable(value) or len(value) != (ncols if side in 'bt' else nrows):
                raise ValueError(f'Need {nmax}-length list of integers for "{side}array" figure with {nmax} {name}, got {side}array={value}.')
            kwargs[side + 'array'] = [*value] # must be listk we test truthiness of contents

    #-------------------------------------------------------------------------#
    # Shared and spanning axes, panel syncing
    #-------------------------------------------------------------------------#
    # Figure out rows and columns "spanned" by each axes in list, for
    # axis sharing and axis label spanning settings
    sharex = int(_default(share, sharex, rc['share']))
    sharey = int(_default(share, sharey, rc['share']))
    if sharex not in range(4) or sharey not in range(4):
        raise ValueError(f'Axis sharing level can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels). Got sharex={sharex} and sharey={sharey}.')
    spanx = _default(span, spanx, 0 if sharex == 0 else None, rc['span'])
    spany = _default(span, spany, 0 if sharey == 0 else None, rc['span'])
    # Get some axes properties, where locations are sorted by axes id.
    # NOTE: These ranges are endpoint exclusive, like a slice object!
    axids = [np.where(array == i) for i in np.sort(np.unique(array)) if i > 0] # 0 stands for empty
    yrange = np.array([[y.min(), y.max()+1] for y,_ in axids]) # range accounting for panels
    xrange = np.array([[x.min(), x.max()+1] for _,x in axids])
    xref = xrange[ref-1,:] # range for reference axes
    yref = yrange[ref-1,:]

    #-------------------------------------------------------------------------#
    # Get basemap.Basemap or cartopy.CRS instances for map, and
    # override aspect ratio.
    #-------------------------------------------------------------------------#
    # NOTE: Cannot have mutable dict as default arg, because it changes the
    # "default" if user calls function more than once! Swap dicts for None.
    basemap = _axes_dict(naxs, basemap, kw=False, default=False)
    proj    = _axes_dict(naxs, _default(projection, proj), kw=False, default='cartesian')
    proj_kw = _axes_dict(naxs, _default(projection_kw, proj_kw, {}), kw=True)
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
            instance, iaspect, kwproj = projs.Proj(name, basemap=basemap[num], **proj_kw[num])
            if num == ref:
                aspect = iaspect
            axes_kw[num].update({'projection':package, 'map_projection':instance})
            axes_kw[num].update(kwproj)

    #-------------------------------------------------------------------------#
    # Figure architecture
    #-------------------------------------------------------------------------#
    # Figure and/or average axes dimensions
    names, values = (), ()
    if journal:
        figsize = _journals(journal) # if user passed width=<string > , will use that journal size
        spec = f'journal={repr(journal)}'
        names = ('axwidth', 'axheight', 'width')
        values = (axwidth, axheight, width)
        width, height = figsize
    elif figsize:
        spec = f'figsize={repr(figsize)}'
        names = ('axwidth', 'axheight', 'width', 'height')
        values = (axwidth, axheight, width, height)
        width, height = figsize
    elif width is not None or height is not None:
        spec = []
        if width is not None:
            spec.append(f'width={repr(width)}')
        if height is not None:
            spec.append(f'height="{repr(height)}"')
        spec = ', '.join(spec)
        names = ('axwidth', 'axheight')
        values = (axwidth, axheight)
    # Raise warning
    for name,value in zip(names,values):
        if value is not None:
            warnings.warn(f'You specified both {spec} and {name}={repr(value)}. Ignoring "{name}".')
    # Standardize input
    width  = units(width)
    height = units(height)
    axwidth  = units(axwidth)
    axheight = units(axheight)

    # Gridspec defaults
    # NOTE: Ratios are scaled to take physical units in _subplots_geometry, so
    # user can manually provide hspace and wspace in physical units.
    wratios = np.atleast_1d(_default(width_ratios, wratios, axwidths, 1))
    hratios = np.atleast_1d(_default(height_ratios, hratios, axheights, 1))
    # Subplot space
    hspace = np.atleast_1d(_default(units(hspace),
        units(rc['subplots.titlespace']) + units(rc['subplots.innerspace']) if sharex == 3
        else units(rc['subplots.xlabspace']) if sharex in (1,2) # space for tick labels and title
        else units(rc['subplots.titlespace']) + units(rc['subplots.xlabspace']
        )))
    wspace = np.atleast_1d(_default(units(wspace),
        units(rc['subplots.innerspace']) if sharey == 3
        else units(rc['subplots.ylabspace']) - units(rc['subplots.titlespace']) if sharey in (1,2) # space for tick labels only
        else units(rc['subplots.ylabspace'])
        ))
    wflush = _default(wflush, flush, False)
    hflush = _default(hflush, flush, False)
    if len(wratios) == 1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios) == 1:
        hratios = np.repeat(hratios, (nrows,))
    if len(wspace) == 1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace) == 1:
        hspace = np.repeat(hspace, (nrows-1,))
    if wflush:
        wspace = wspace*0.0
    if hflush:
        hspace = hspace*0.0
    # Border space
    left   = units(_default(left,   rc['subplots.ylabspace']))
    bottom = units(_default(bottom, rc['subplots.xlabspace']))
    right  = units(_default(right,  rc['subplots.nolabspace']))
    top    = units(_default(top,    rc['subplots.titlespace']))

    # Parse arguments, fix dimensions in light of desired aspect ratio
    figsize, gridspec_kw, subplots_kw = _subplots_geometry(
        nrows=nrows, ncols=ncols, aspect=aspect, xref=xref, yref=yref,
        left=left, right=right, bottom=bottom, top=top,
        width=width, height=height, axwidth=axwidth, axheight=axheight,
        wratios=wratios, hratios=hratios, wspace=wspace, hspace=hspace,
        wpanels=[0]*ncols, hpanels=[0]*nrows, # computed automatically later
        **kwargs)
    # Apply settings and add attributes
    gridspec = FlexibleGridSpec(**gridspec_kw)
    # Create blank figure
    # WARNING: Important to pass main_gridspec and supblots_kw to figure right
    # away because some backends will call fig.draw immediately on instantiation.
    fig = plt.figure(FigureClass=Figure, tight=tight, figsize=figsize, ref=ref,
        tightborders=tightborders, tightsubplots=tightsubplots, tightpanels=tightpanels,
        borderpad=borderpad, subplotpad=subplotpad, panelpad=panelpad,
        flush=flush, wflush=wflush, hflush=hflush,
        autoformat=autoformat,
        main_gridspec=gridspec,
        subplots_kw=subplots_kw,
        )
    fig._locked = False

    #-------------------------------------------------------------------------#
    # Draw on figure
    #-------------------------------------------------------------------------#
    # Draw main subplots
    axs = naxs*[None] # list of axes
    for idx in range(naxs):
        # Get figure gridspec ranges
        num = idx + 1
        x0, x1 = xrange[idx,0]+1, xrange[idx,1]+1
        y0, y1 = yrange[idx,0]+1, yrange[idx,1]+1
        # Get gridspec for ensemble of axes and panels
        # TODO: What happens when figure size changes? Are these fixed?
        subplotspec = gridspec[y0:y1, x0:x1]
        figwidth, figheight = fig.get_size_inches()
        bbox = subplotspec.get_position(fig) # valid since axes not drawn yet
        width = abs(bbox.width)*figwidth
        height = abs(bbox.height)*figheight
        igridspec = FlexibleGridSpecFromSubplotSpec(
                subplot_spec=subplotspec,
                nrows=3, ncols=3,
                wspace=[0,0], hspace=[0,0],
                width_ratios=[0,width,0], height_ratios=[0,height,0],
                )
        # Add subplot
        axs[idx] = fig.add_subplot(igridspec[1,1],
            number=num, spanx=spanx, spany=spany,
            sharex_level=sharex, sharey_level=sharey,
            main_gridspec=igridspec,
            **axes_kw[num])
    # Set up shared axes
    fig._main_axes = axs
    if sharex:
        groups = fig._get_shared_axes('x') # groups sharing unique x-extents in the gridspec
        for group in groups:
            parent = group[0]
            for child in group[1:]:
                child._sharex_setup(parent, sharex)
    if sharey:
        groups = fig._get_shared_axes('y')
        for group in groups:
            parent = group[0]
            for child in group[1:]:
                child._sharey_setup(parent, sharey)
    # Draw axes panels after all main subplots are drawn
    # NOTE: This *must* come after shared axes are set up! Otherwise tight
    # layout scaling are wrong.
    for idx in range(naxs):
        num = idx + 1
        ax = axs[idx]
        kw = axpanels_kw[num]
        sides = kw.pop('sides')
        for side in sides:
            ax.panel_axes(side, order=order, **kw)

    # Create outer panels
    for side in 'blrt':
        # Draw panel from gridspec
        array = subplots_kw[side + 'array']
        if array is None:
            continue
        name = _side_translate[side]
        paxs = []
        for num in np.unique(array).flat:
            # Get subspec, and settings
            # TODO: Allow for different width ratios and stuff
            if num == 0:
                continue
            flush = kwargs[side + 'flush']
            idx, = np.where(array == num)
            idx = slice(min(idx)+1, max(idx)+2)
            if side in 'lr':
                hspace, hratios = [], 1
                wspace, wratios = kwargs[side + 'sep'], kwargs[side + 'width']
                nrows, ncols = 1, len(wratios)
                if side == 'r':
                    subspec = gridspec[idx,-1] # ignore corners
                else:
                    subspec = gridspec[idx,0]
            else:
                wspace, wratios = [], 1
                hspace, hratios = kwargs[side + 'sep'], kwargs[side + 'width']
                nrows, ncols = len(hratios), 1
                if side == 'b':
                    subspec = gridspec[-1,idx]
                else:
                    subspec = gridspec[0,idx]
            # Make gridspec for containing the "stack" of panels
            ipaxs = []
            igridspec = FlexibleGridSpecFromSubplotSpec(subplot_spec=subspec,
                nrows=nrows, ncols=ncols,
                wspace=wspace, hspace=hspace,
                width_ratios=wratios, height_ratios=hratios,
                )
            for i in range(max((nrows,ncols))):
                ipax = fig.add_subplot(igridspec[i], projection='panel', side=name, flush=flush)
                ipax._main_gridspec = fig._main_gridspec
                ipax._stack_gridspec = igridspec
                ipaxs += [ipax]
            paxs += [ipaxs]
        # Sort panel axes into row-major or column-major order
        if (side in 'bt' and order == 'C') or (side in 'lr' and order != 'C'):
            paxs = [*zip(*paxs)]
        # Store in axes_grid with support for 2D indexing
        n = len(paxs[0])
        paxs = [ax for ipaxs in paxs for ax in ipaxs] # unfurl
        paxs = axes_grid(paxs, n=n, order=order)
        setattr(fig, '_' + name + 'panel', paxs)

    # Return figure and axes
    n = (ncols if order == 'C' else nrows)
    fig._locked = True
    return fig, axes_grid(axs, n=n, order=order)

