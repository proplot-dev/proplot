#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures and axes.
The `subplots` function is all you'll need to directly use here.
It returns a `Figure` instance and an `axes_grid` container of
`~proplot.axes.BaseAxes` axes, whose positions are controlled by the
`FlexibleGridSpec` class.

.. raw:: html

   <h1>Developer notes</h1>

Many of ProPlot's more complex features are stuffed right into the `subplots`
function. But why didn't we separate these features into their own functions --
for example, ``ax.add_panel('right')`` or ``ax.set_projection('proj')``?

ProPlot requires a *static* grid of subplots and panels before anything is plotted,
so that we can exert more control on the layout and make things look "nice"
without any manual tweaking on the part of the user. `subplots` automatically
aligns axes and figure panels in the subplot grid and controls subplot aspect
ratios. `~Figure.smart_tight_layout` automatically adjusts aspect ratios for map
projection and imshow plots, and controls the amount of whitespace between subplot
content, panel content, and the figure edge.
"""
import os
import re
import numpy as np
from numbers import Number
import functools
import warnings
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.figure as mfigure
import matplotlib.transforms as mtransforms
import matplotlib.gridspec as mgridspec
from .rctools import rc
from .utils import _default, units
from . import projs, axes
__all__ = [
    'axes_grid', 'close', 'show', 'subplots', 'Figure',
    'FlexibleGridSpec', 'FlexibleGridSpecBase', 'FlexibleGridSpecFromSubplotSpec',
    ]

# Aliases for panel names
_panel_aliases = {
    'bpanel':         'bottompanel',
    'rpanel':         'rightpanel',
    'lpanel':         'leftpanel',
    'bcolorbar':      'bottompanel',
    'rcolorbar':      'rightpanel',
    'lcolorbar':      'leftpanel',
    'blegend':        'bottompanel',
    'rlegend':        'rightpanel',
    'llegend':        'leftpanel',
    'bottomcolorbar': 'bottompanel',
    'rightcolorbar':  'rightpanel',
    'leftcolorbar':   'leftpanel',
    'bottomlegend':   'bottompanel',
    'rightlegend':    'rightpanel',
    'leftlegend':     'leftpanel',
    }

#------------------------------------------------------------------------------#
# Miscellaneous stuff
#------------------------------------------------------------------------------#
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
    stacked axes panels. See the `~axes_grid.__getattr__` and
    `~axes_grid.__getitem__` methods for details."""
    def __init__(self, list_, n=1, order='C'):
        # Add special attributes that support 2D grids of axes
        # NOTE: The input list is always a vector *already unfurled* in row-major
        # or column-major order, and 'n' is the fastest-moving dimension size, i.e.
        # ncols for order=='C' and nrows for order=='F'.
        self._n = n # ncols or nrows
        self._order = order # order
        super().__init__(list_)

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
        if isinstance(key, tuple) and len(key)==1:
            key = key[0]
        if not isinstance(key, tuple): # do not expand single slice to list of integers or we get recursion! len() operator uses __getitem__!
            axlist = isinstance(key, slice)
            objs = list.__getitem__(self, key)
        elif len(key)==2:
            axlist = any(isinstance(ikey, slice) for ikey in key)
            # Expand keys
            # Testing shows that this order stuff works.
            keys = []
            order = self._order
            for i,ikey in enumerate(key):
                if isinstance(ikey, slice):
                    start, stop, step = ikey.start, ikey.stop, ikey.step
                    if start is None:
                        start = 0
                    if step is None:
                        step = 1
                    if stop is None:
                        if (i==1 and order=='C') or (i==0 and order!='C'):
                            stop = self._n
                        else:
                            stop = len(self)//self._n
                    ikeys = [*range(start, stop, step)]
                else:
                    ikeys = [ikey]
                keys.append(ikeys)
            # Get index pairs and get objects
            # Note that in double for loop, right loop varies fastest, so
            # e.g. axs[:,:] delvers (0,0), (0,1), ..., (0,N), (1,0), ...
            # Remember for order=='F', axes_grid was sent a list unfurled in
            # column-major order, so we replicate row-major indexing syntax by
            # reversing the order of the keys.
            objs = []
            if self._order=='C':
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
        on multiple axes at once! If the attribute is *not callable*, returns
        an `axes_grid` of identically named attributes for every object in the
        list.

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
        # Objects
        if not any(callable(_) for _ in attrs):
            return axes_grid(attrs) # or just attrs?
        # Methods
        if all(callable(_) for _ in attrs):
            @functools.wraps(attrs[0])
            def axes_grid_iterator(*args, **kwargs):
                ret = []
                for func in attrs:
                    ret.append(func(*args, **kwargs))
                if len(ret)==1:
                    return ret[0]
                elif all(res is None for res in ret):
                    return None
                else:
                    return (*ret,) # expand to tuple
            return axes_grid_iterator
        # Mixed
        raise AttributeError(f'Found mixed types for attribute "{attr}".')

#------------------------------------------------------------------------------#
# Gridspec classes
#------------------------------------------------------------------------------#
def _adjust(n):
    """Account for negative indices."""
    if n<0:
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
        # gridspec (see 'update')
        self._nrows_visible = nrows
        self._ncols_visible = ncols
        self._nrows = nrows*2-1
        self._ncols = ncols*2-1
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
        nrows, ncols = self._nrows, self._ncols
        nrows_visible, ncols_visible = self._nrows_visible, self._ncols_visible
        if not isinstance(key, tuple): # usage gs[1,2]
            num1, num2 = _normalize(key, nrows_visible * ncols_visible)
        else:
            if len(key)==2:
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
        nrows = self._nrows_visible
        ncols = self._ncols_visible
        hratios = np.atleast_1d(_default(height_ratios, 1))
        wratios = np.atleast_1d(_default(width_ratios,  1))
        hspace = np.atleast_1d(_default(hspace, np.mean(hratios)*0.10)) # this is relative to axes
        wspace = np.atleast_1d(_default(wspace, np.mean(wratios)*0.10))
        if len(wspace)==1:
            wspace = np.repeat(wspace, (ncols-1,)) # note: may be length 0
        if len(hspace)==1:
            hspace = np.repeat(hspace, (nrows-1,))
        if len(wratios)==1:
            wratios = np.repeat(wratios, (ncols,))
        if len(hratios)==1:
            hratios = np.repeat(hratios, (nrows,))

        # Verify input ratios and spacings
        # Translate height/width spacings, implement as extra columns/rows
        if len(hratios) != nrows:
            raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')
        if len(wratios) != ncols:
            raise ValueError(f'Got {ncols} columns, but {len(wratios)} wratios.')
        if ncols>1 and len(wspace) != ncols-1:
            raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(wspace)}.')
        if nrows>1 and len(hspace) != nrows-1:
            raise ValueError(f'Require {nrows-1} height spacings for {nrows} rows, got {len(hspace)}.')

        # Assign spacing as ratios
        # Also return extra kwargs, will be passed to superclass initializers
        # or superclass update method.
        wratios_final = [None]*self._ncols
        wratios_final[::2] = [*wratios]
        if self._ncols>1:
            wratios_final[1::2] = [*wspace]
        hratios_final = [None]*self._nrows
        hratios_final[::2] = [*hratios]
        if self._nrows>1:
            hratios_final[1::2] = [*hspace]
        return wratios_final, hratios_final, kwargs # bring extra kwargs back

    def update(self, **gridspec_kw):
        """Update the width, height ratios and spacing for subplot columns, rows."""
        wratios, hratios, edges_kw = self._spaces_as_ratios(**gridspec_kw)
        self.set_width_ratios(wratios)
        self.set_height_ratios(hratios)
        edges_kw.pop('nrows', None) # cannot be modified
        edges_kw.pop('ncols', None)
        super().update(**edges_kw) # remaining kwargs should just be left, right, top, bottom

class FlexibleGridSpec(FlexibleGridSpecBase, mgridspec.GridSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpec`."""
    pass

class FlexibleGridSpecFromSubplotSpec(FlexibleGridSpecBase, mgridspec.GridSpecFromSubplotSpec):
    """Mixes `FlexibleGridSpecBase` with `~matplotlib.gridspec.GridSpecFromSubplotSpec`."""
    pass

#------------------------------------------------------------------------------#
# Figure class
#------------------------------------------------------------------------------#
def _iter_twins(ax):
    """Iterates over twin axes."""
    # TODO: should this include inset axes? or should we ignore inset
    # axes when doing tight layouts?
    axs = []
    if not ax:
        return axs
    if not ax.get_visible():
        return axs
    for iax in (ax, ax._altx_child, ax._alty_child):
        if not iax:
            continue
        if not iax.get_visible():
            continue
        axs.append(iax)
    return axs

def _iter_children(ax):
    """Iterates over an axes and its children, including panels, twin
    axes, and panel twin axes"""
    axs = []
    if not ax:
        return axs
    if not ax.get_visible():
        return axs
    iaxs = (ax, *ax.leftpanel, *ax.bottompanel, *ax.rightpanel, *ax.toppanel)
    for iax in iaxs:
        axs.extend(_iter_twins(iax))
    return axs

def _intervalx(ax):
    """Given an axes and a bounding box, pads the intervalx according to the
    matplotlib "tight layout" error associated with invisible y ticks."""
    bbox = ax._tight_bbox
    if not isinstance(ax, axes.CartesianAxes):
        return (bbox.xmin, bbox.xmax)
    xerr = ax._ytick_pad_error # error in x-direction, due to y ticks
    return (bbox.xmin + xerr[0], bbox.xmax - sum(xerr))

def _intervaly(ax):
    """Given an axes and a bounding box, pads the intervaly according to the
    matplotlib "tight layout" error associated with invisible x ticks."""
    bbox = ax._tight_bbox
    if not isinstance(ax, axes.CartesianAxes):
        return (bbox.ymin, bbox.ymax)
    yerr = ax._xtick_pad_error # error in y-direction, due to x ticks
    return (bbox.ymin + yerr[0], bbox.ymax - sum(yerr))

def _ax_span(ax, renderer, children=True):
    """Get span, accounting for panels, shared axes, and whether axes has
    been replaced by colorbar in same location."""
    # Return arrays
    axs = _iter_children(ax)
    axs = [ax for ax in axs if ax._tight_bbox is not None]
    xs = np.array([_intervalx(ax) for ax in axs])
    ys = np.array([_intervaly(ax) for ax in axs])
    xspan = [xs[:,0].min(), xs[:,1].max()]
    yspan = [ys[:,0].min(), ys[:,1].max()]
    return xspan, yspan

def _ax_props(axs, renderer):
    """If this is a panel axes, check if user generated a colorbar axes,
    which stores all the artists/is what we really want to get a
    tight bounding box for."""
    axs = [ax for ax in axs if (ax and ax.get_visible())]
    if not axs:
        return (np.empty((0,2)),)*4
    spans = [_ax_span(ax, renderer) for ax in axs]
    xspans = np.array([span[0] for span in spans])
    yspans = np.array([span[1] for span in spans])
    xrange = np.array([(ax._colorbar_parent or ax)._xrange for ax in axs])
    yrange = np.array([(ax._colorbar_parent or ax)._yrange for ax in axs])
    return xspans, yspans, xrange, yrange

class Figure(mfigure.Figure):
    def __init__(self,
            tight=True, tightborders=None, tightsubplots=None, tightpanels=None,
            flush=False, wflush=None, hflush=None,
            borderpad=None, subplotpad=None, panelpad=None,
            autoformat=True,
            **kwargs):
        """
        The `~matplotlib.figure.Figure` instance returned by `subplots`. At
        draw-time, this class conveniently adjusts subplot positioning and the
        outer figure bounding box to accomodate text and plotted content,
        without messing up subplot aspect ratios and panel widths. See
        `~Figure.smart_tight_layout` for details.

        Parameters
        ----------
        tight : bool, optional
            Default setting for `tightborders`, `tightsubplots`, and `tightpanels`
            keyword args. Defaults to ``rc['tight']``.
        tightborders : bool, optional
            Whether to draw a tight bounding box around the whole figure.
            If ``None``, takes the value of `tight`.
        tightsubplots : bool, optional
            Whether to automatically space out subplots to prevent overlapping
            axis tick labels, etc. If ``None``, takes the value of `tight`.
        tightpanels : bool, optional
            Whether to automatically space between subplots and their panels
            to prevent overlap. If ``None``, takes the value of `tight`.
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

        Other parameters
        ----------------
        **kwargs
            Passed to `matplotlib.figure.Figure`.
        """
        # Tight toggling
        tight = _default(tight, rc['tight'])
        self._smart_tight_outer   = _default(tightborders, tight)
        self._smart_tight_subplot = _default(tightsubplots, tight)
        self._smart_tight_panel   = _default(tightpanels, tight)
        self._smart_tight_init = True # is figure in its initial state?
        # Padding and "flush" args
        self._extra_pad = 0 # sometimes matplotlib fails, cuts off super title! will add to this
        self._smart_borderpad = units(_default(borderpad, rc['subplots.borderpad']))
        self._smart_subplotpad  = units(_default(subplotpad,  rc['subplots.subplotpad']))
        self._smart_panelpad = units(_default(panelpad, rc['subplots.panelpad']))
        self._subplot_wflush = _default(wflush, flush)
        self._subplot_hflush = _default(hflush, flush)
        # Gridspec information, filled in by subplots()
        self._ref_num = 1
        self._main_axes = []  # list of 'main' axes (i.e. not insets or panels)
        self._spanning_axes = [] # add axis instances to this, and label position will be updated
        self._subplots_kw = None # extra special settings
        self._main_gridspec = None # gridspec encompassing drawing area
        # Figure-wide settings
        self._locked = True
        self._autoformat = autoformat
        # Panels, initiate as empty
        self.leftpanel   = axes.EmptyPanel()
        self.bottompanel = axes.EmptyPanel()
        self.rightpanel  = axes.EmptyPanel()
        self.toppanel    = axes.EmptyPanel()
        # Initialize
        super().__init__(**kwargs)
        self.suptitle('') # this adds _suptitle position
        self._suptitle_transform = None

    def __getattribute__(self, attr, *args):
        """Enables the attribute aliases ``bpanel`` for ``bottompanel``,
        ``tpanel`` for ``toppanel``, ``lpanel`` for ``leftpanel``, and
        ``rpanel`` for ``rightpanel``. Issues warning for new users that
        try to access the `~matplotlib.figure.Figure.add_subplot` and
        `~matplotlib.figure.Figure.colorbar` functions."""
        # Just test for add_subplot, do not care if user provides figure
        # colorbar method with an axes manually.
        attr = _panel_aliases.get(attr, attr)
        if attr=='add_subplot' and self._locked:
            warnings.warn('Using "add_subplot" or "colorbar" with ProPlot figures may result in unexpected behavior. '
                'Please use the subplots() command to create your figure, subplots, and panels all at once.')
        return super().__getattribute__(attr, *args)

    def _suptitle_setup(self, title, **kwargs):
        """Assign figure "super title"."""
        if title is not None and self._suptitle.get_text()!=title:
            self._suptitle.set_text(title)
        if kwargs:
            self._suptitle.update(kwargs)

    def _labels_setup(self, ax, labels, rows=True, **kwargs):
        """Assigns row labels and column labels, updates label settings."""
        attr = ('rowlabel' if rows else 'collabel')
        range_along = ('_yrange' if rows else '_xrange') # along labels
        range_across = ('_xrange' if rows else '_yrange') # across from labels
        if not self._main_axes:
            axs = [ax]
        else:
            axs = [ax for ax in self._main_axes if getattr(ax, range_across)[0]==0] # axes on the edge
            axs = [ax for _,ax in sorted(zip([getattr(ax, range_along)[0] for ax in axs], axs))] # order by yrange
        if labels is None or isinstance(labels, str): # common during testing
            labels = [labels]*len(axs)
        if len(labels)!=len(axs):
            raise ValueError(f'Got {len(labels)} {attr}s, but there are {len(axs)} {attr.split("label")[0]}s.')
        for ax,label in zip(axs,labels):
            obj = getattr(ax, attr)
            if label is not None and obj.get_text()!=label:
                obj.set_text(label)
            if kwargs:
                obj.update(kwargs)

    def _tight_bboxs(self, renderer=None):
        """Sets the ``_tight_bbox`` attribute on axes, so that we don't have
        to call the tight bbox algorithm multiple times."""
        for ax in (*self._main_axes, *self.leftpanel, *self.bottompanel, *self.rightpanel):
            for iax in _iter_children(ax):
                bbox = (iax._colorbar_child or iax).get_tightbbox(renderer)
                iax._tight_bbox = bbox

    def _post_aspect_fix(self):
        """Adjust average aspect ratio used for gridspec calculations."""
        # This will fix grids with identically fixed aspect ratios, i.e.
        # identically zoomed-in cartopy projections or imshow images. For more
        # complex grids with varying fixed aspect ratios, extra whitespace will
        # be inevitable unless user fiddles with hratios and wratios.
        if not self._main_axes:
            return
        ax = self._main_axes[self._ref_num-1]
        aspect = None
        subplots_kw = self._subplots_kw
        if isinstance(ax, axes.CartesianAxes):
            aspect = ax._aspect_equal
            ax._aspect_equal = None
        elif isinstance(ax, axes.CartopyProjectionAxes):
            bbox = ax.background_patch._path.get_extents()
            aspect = abs(bbox.width)/abs(bbox.height)
            if aspect==subplots_kw['aspect']:
                aspect = None
        if aspect is not None:
            subplots_kw['aspect'] = aspect
            figsize, gridspec_kw, _ = _subplots_kwargs(**subplots_kw)
            self._main_gridspec.update(**gridspec_kw)
            self.set_size_inches(figsize)

    def _axis_label_update(self, axis, span=False, **kwargs):
        """Get axis label for axes with axis sharing or spanning enabled.
        When `span` is False, we add labels to all axes. When `span` is
        True, we filter to one axes."""
        # Get the axes
        # TODO: Account for axis with different size tick labels?
        # NOTE: Only do this *after* tight bounding boxes have been drawn!
        name = axis.axis_name
        if name not in 'xy': # i.e. theta or radius
            return
        base = axis.axes
        for _ in range(2): # try 2 levels down, should be sufficient
            base = getattr(base, '_share' + name, None) or base
        span_on = getattr(base, '_span'  + name)
        if not span_on:
            axes = [getattr(base, name + 'axis')]
        else:
            if name=='x':
                axs = [ax for ax in self._main_axes if ax._yrange[1]==base._yrange[1]]
            else:
                axs = [ax for ax in self._main_axes if ax._xrange[0]==base._xrange[0]]
            axes = [getattr(getattr(ax, '_share' + name) or ax, name + 'axis') for ax in axs]

        # Apply settings and store list of spanning axes
        for axis in axes:
            if span_on and axis not in self._spanning_axes:
                self._spanning_axes.append(axis)
            if axis.get_label_position() == 'top':
                kwargs['va'] = 'bottom' # baseline gets cramped if no ticklabels present
            axis.label.update({'visible':True, **kwargs})

        # If requested, turn on spanning
        # TODO: Use tightbbox for this?
        if span and span_on:
            # Use halfway-point axis for spanning label, because if we have an
            # odd number of columns or rows, that transform gives us perfect offset
            ranges = np.array([getattr(ax, '_' + name + 'range') for ax in axs])
            half = (np.argmin(ranges[:,0]) + np.argmax(ranges[:,0]))//2
            saxis = axes[half]
            for axis in axes:
                if saxis is not axis:
                    axis.label.update({'visible':False})
            # Reposition to "span" other axes
            idx = slice(ranges.min(), ranges.max() + 1)
            if name=='x': # span columns
                subspec = self._main_gridspec[0,idx]
            else: # spans rows
                subspec = self._main_gridspec[idx,0]
            bbox = subspec.get_position(self) # in figure-relative coordinates
            x0, y0, width, height = bbox.bounds
            if name=='x':
                transform = mtransforms.blended_transform_factory(self.transFigure, mtransforms.IdentityTransform())
                position = (x0 + width/2, 1)
            else:
                transform = mtransforms.blended_transform_factory(mtransforms.IdentityTransform(), self.transFigure)
                position = (1, y0 + height/2)
            saxis.label.update({'position':position, 'transform':transform})

    def _post_process_misc(self):
        """Does various post-processing steps required due to user actions
        after the figure was created."""
        # Apply various post processing on per-axes basis
        axs = (iax for ax in self._main_axes for iax in _iter_children(ax))
        for ax in axs:
            if not ax:
                continue
            elif not ax.get_visible():
                continue
            # Lock dual axes limits
            for xy in 'xy':
                scale = getattr(ax, f'_dual{xy}_scale')
                if not scale:
                    continue
                xyscale = getattr(ax, f'get_{xy}scale')()
                olim = getattr(ax, f'get_{xy}lim')()
                twin = getattr(ax, f'_alt{xy}_child')
                transform = getattr(twin, f'{xy}axis')._scale.get_transform()
                if re.match('^log', xyscale) and any(np.array(olim)<=0):
                    raise ValueError('Axis limits go negative, and "alternate units" axis uses log scale.')
                elif re.match('^inverse', xyscale) and any(np.array(olim)<=0):
                    raise ValueError('Axis limits cross zero, and "alternate units" axis uses inverse scale.')
                # Apply new axis lim
                lim = transform.inverted().transform(np.array(olim))
                if np.sign(np.diff(olim)) != np.sign(np.diff(lim)): # the transform flipped it, so when we try to set limits, will get flipped again!
                    lim = lim[::-1]
                getattr(twin, f'set_{xy}lim')(scale[0] + scale[1]*lim)
            # Axis rotation, if user drew anything that triggered datetime x-axis
            # WARNING: Rotation is done *before* horizontal/vertical alignment, and
            # cannot change alignment with set_tick_params. Need to iterate through
            # text objects. Also do not use fig.autofmt_date, it calls subplots_adjust
            if not ax._xrotated and isinstance(ax.xaxis.converter, mdates.DateConverter):
                rotation = rc['axes.formatter.timerotation']
                kw = {'rotation':rotation}
                if rotation not in (0,90,-90):
                    kw['ha'] = ('right' if rotation>0 else 'left')
                for label in ax.xaxis.get_ticklabels():
                    label.update(kw)
            # Automatic labels and colorbars for plot
            # NOTE: The legend wrapper supports multiple legends in one axes
            # by adding legend artists manually.
            for loc,handles in ax._auto_colorbar.items():
                ax.colorbar(handles, **ax._auto_colorbar_kw[loc])
            for loc,handles in ax._auto_legend.items():
                ax.legend(handles, **ax._auto_legend_kw[loc]) # deletes other ones!

    def _post_text_align(self, renderer):
        """Adjusts position of row titles and figure super title."""
        # Adjust row labels as axis tick labels are generated and y-axis
        # labels is generated.
        # WARNING: Must call get_tightbbox manually, since the _tight_bbox
        # attribute *includes* labels for appropriate tight layout, but we need
        # to *exlude* labels before positioning them by making them invisible!
        w, h = self.get_size_inches()
        axs = [ax for ax in self._main_axes if ax._xrange[0]==0] # order by xrange
        for ax in axs:
            label = ax.rowlabel
            label.set_visible(False) # make temporarily invisible, so tightbbox does not include existing label!
            if not label.get_text().strip():
                continue
            # Iterate panels
            ixs = []
            if ax.leftpanel:
                iaxs = ax.leftpanel
            else:
                iaxs = (ax, *ax.toppanel, *ax.bottompanel)
            for iax in iaxs:
                for jax in _iter_twins(iax):
                    # Box for column label
                    # bbox = jax._tight_bbox
                    bbox = jax.get_tightbbox(renderer)
                    x, _ = self.transFigure.inverted().transform((bbox.xmin, 0))
                    ixs.append(x)
            # Verify, be careful!
            label.set_visible(True)
            if not ixs:
                warnings.warn('Axes on left row are invisible. Cannot determine rowtitle position.')
                continue
            # Update position
            # bbox = ax._tight_bbox
            label.set_visible(True)
            x = min(ixs)
            transform = mtransforms.blended_transform_factory(self.transFigure, ax.transAxes)
            x = x - (0.6*label.get_fontsize()/72)/w # see: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
            label.update({'x':x, 'y':0.5, 'ha':'right', 'va':'center', 'transform':transform})

        # Adjust col labels -- this is much simpler
        ys = []
        suptitle = self._suptitle
        suptitle.set_visible(False)
        axs = [ax for ax in self._main_axes if ax._yrange[0]==0] # order by xrange
        for ax in axs:
            label = ax.collabel
            label.set_visible(False)
            if not label.get_text().strip() and not suptitle.get_text().strip():
                continue
            # Get maximum y-position among all children, necessary because
            # e.g. twin x-axes with their own labels are common
            iys = []
            if ax.toppanel and ax.toppanel.get_visible():
                iaxs = ax.toppanel
            else:
                iaxs = (ax, *ax.leftpanel, *ax.rightpanel)
            for iax in iaxs:
                for jax in _iter_twins(iax):
                    # Box for column label
                    # bbox = jax._tight_bbox
                    bbox = jax.get_tightbbox(renderer)
                    _, y = self.transFigure.inverted().transform((0, bbox.ymax))
                    iys.append(y)
            # Update column label position
            label.set_visible(True)
            if not iys:
                warnings.warn('Axes on top row is invisible. Cannot determine coltitle position.')
                continue
            if label.get_text().strip():
                y = max(iys)
                transform = mtransforms.blended_transform_factory(ax.transAxes, self.transFigure)
                y = y + (0.3*label.get_fontsize()/72)/h # see: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
                label.update({'x':0.5, 'y':y, 'ha':'center', 'va':'bottom', 'transform':transform})
            # Box for super title
            # This time account for column labels of course!
            if suptitle.get_text().strip():
                if label.get_text().strip(): # by construction it is above everything else!
                    # bbox = ax._tight_bbox
                    bbox = ax.get_tightbbox(renderer)
                    _, y = self.transFigure.inverted().transform((0, bbox.ymax))
                else:
                    y = max(iys)
                ys.append(y)

        # Super title position
        # Get x position as center between left edge of leftmost main axes
        # and right edge of rightmost main axes
        suptitle.set_visible(True)
        if suptitle.get_text().strip():
            if not ys:
                warnings.warn('All axes on top row are invisible. Cannot determine suptitle position.')
            else:
                kw = self._subplots_kw
                left = kw['left']
                right = kw['right']
                if self.leftpanel:
                    left += (kw['lwidth'] + kw['lspace'])
                if self.rightpanel:
                    right += (kw['rwidth'] + kw['rspace'])
                x = left/w + 0.5*(w - left - right)/w
                y = max(ys) + (0.3*suptitle.get_fontsize()/72)/h
                suptitle.update({'x':x, 'y':y, 'ha':'center', 'va':'bottom', 'transform':self.transFigure})

    def _post_process(self, renderer=None):
        """Performs various post-procesing tasks."""
        # Get renderer
        if renderer is None:
            renderer = self.canvas.get_renderer()
        # Post-plotting stuff
        if self._smart_tight_init:
            self._post_aspect_fix()
            self._post_process_misc()
        # Row, column, figure, spanning labels
        # WARNING: draw() is called *more than once* and title positions are
        # appropriately offset only during the *later* calls! Must run each time.
        for axis in self._spanning_axes: # turn spanning off so rowlabels adjust properly
            self._axis_label_update(axis, span=False)
        self._post_text_align(renderer) # just applies the spacing
        # Tight layout
        # WARNING: For now just call this once, get bugs if call it every
        # time draw() is called, but also means bbox accounting for titles
        # and labels may be slightly off
        if not self._smart_tight_init or not (self._smart_tight_outer or self._smart_tight_subplot or self._smart_tight_panel):
            pass
        elif any(ax._gridliner_on for ax in self._main_axes):
            warnings.warn('Tight subplots do not work with cartopy gridline labels. Use tight=False in your call to subplots().')
        else:
            self.smart_tight_layout(renderer)
        for axis in self._spanning_axes: # turn spanning back on
            self._axis_label_update(axis, span=True)
        # Set flag
        self._smart_tight_init = False

    def _panel_tight_layout(self, side, paxs, renderer, figure=False):
        """From a list of panels axes int he same row or column, figure out
        the necessary new 'spacing' to apply to the inner gridspec object.
        For axes panels, this function just modifies mutable width and height
        ratio lists in place, and returns the 'wpanels' or 'hpanels' arguments.
        For figure panels, it returns the 'lsep', 'rsep', or 'bsep' argument
        needed for _subplots_kwargs."""
        # Initial stuff
        # Check that none of input panels are EmptyPanel, although they can
        # certainly be invisible/allocated to make spacing even.
        pad = self._smart_panelpad
        if not paxs:
            return 0
        elif all(not ipanel for pax in paxs for ipanel in pax): # exists, but may be invisible
            return 0
        elif not all(ipanel for pax in paxs for ipanel in pax):
            raise ValueError('Either all or no axes in this row should have panels.')
        elif len({len(pax) for pax in paxs})>1:
            raise ValueError(f'Different numbers of stacked panels in same row or column: {[len(pax) for pax in paxs]}.')
        # Get list of panels in each panel stack that are adjacent to the
        # main subplot area.
        if side in 'br':
            paxs_main = [pax[0] for pax in paxs]
        else:
            paxs_main = [pax[-1] for pax in paxs]
        # Get full width and height ratios for the panel and outer gridspecs
        # NOTE: We edit the *ratio lists* themselves. The getters do not
        # return copies, but the mutable underlying objects!
        gridspecs = [pax._panels_main_gridspec for pax in paxs_main]
        pgridspecs = [pax._panels_stack_gridspec for pax in paxs_main]
        if side in 'lr':
            ratios = [gs.get_width_ratios() for gs in gridspecs]
            pratios = [gs.get_width_ratios() for gs in pgridspecs]
        else:
            ratios = [gs.get_height_ratios() for gs in gridspecs]
            pratios = [gs.get_height_ratios() for gs in pgridspecs]

        # Get necessary space between stacked panels. We iterate through pairs
        # of panels and get *minimum* actual space between adjacent panels.
        seps = []
        for pax in paxs:
            # Ratios
            isep = []
            for i in range(len(pax)-1): # empty if there is just one panel in the stack
                ipaxs = pax[i:i+2] # lists are left-to-right, top-to-bottom
                if any(ipax._tight_bbox is None for ipax in ipaxs):
                    isep.append(None)
                    continue
                # Get intervals and spacing. Options are:
                # 1) Bottom of top panel minus top of bottom panel
                # 2) Left of right panel minus right of left panel
                # NOTE: Order of stacks is always left-right and top-bottom
                if side in 'lr':
                    ispans = [_intervalx(pax) for pax in ipaxs]
                    isep.append((ispans[1][0] - ispans[0][1])/self.dpi) # bottom of top one minus top of bottom one
                else: # 'tb'
                    ispans = [_intervaly(pax) for pax in ipaxs]
                    isep.append((ispans[0][0] - ispans[1][1])/self.dpi) # bottom of top one minus top of bottom one
            # List of actual spacing between panels in each panel stack
            seps.append(isep)

        # Setup
        sep = []
        flush = [pax._flush for pax in paxs_main]
        if not all(flush) and any(flush):
            warnings.warn(f'"flush" setting varies between panels in the same row or column. Using "flush" from the first axes in the row or column: {flush[0]}.')
        flush = flush[0]
        # Apply updated *space* to the underlying panel gridspec and main
        # axes gridspec width and height ratios
        # NOTE: The existing seps in ipratios will be *identical* along
        # the row or column, so we use a sample.
        seps = [[*_] for _ in zip(*seps)] # seps grouped by like row or column
        for idx,isep in enumerate(seps):
            idx = 1 + idx*2 # index of *spaces* in ratios list
            isep = [i for i in isep if i is not None] # if panel invisible
            if not isep:
                warnings.warn('All panels in this row or column at some stack level are invisible. That is really weird.')
                sep.append(0)
                continue
            if flush: # always want panels touching, in spite of ticks, etc.
                isep = 0
            else:
                isep = max([0, pratios[0][idx] - min(isep) + pad])
            sep.append(isep)
        # Adjust ratios for panel gridspec
        for iratios,ipratios in zip(ratios,pratios):
            for idx,isep in enumerate(sep):
                idx = 1 + idx*2 # index of *spaces* in ratios list
                ipratios[idx] = isep
                if figure:
                    continue
                idx_stack = (-1 if side in 'br' else 0)
                iratios[idx_stack] = sum(ipratios)
        # Just need the 'sep' argument for outer panels. We have already adjusted
        # 'lspace', etc. arguments by modifying the wratio and hratio lists.
        if figure:
            return sep

        # Get space between main subplot and adjacent panels. Iterate through
        # panels, get maximum necessary spacing. Assign equal spacing to all
        # so they are still aligned
        spaces = []
        for pax in paxs_main:
            if pax._tight_bbox is None:
                continue
            if side in 'lr':
                pspan = _intervalx(pax)
                span = _intervalx(pax._parent)
            else:
                pspan = _intervaly(pax)
                span = _intervaly(pax._parent)
            if side in 'tr':
                spaces.append((pspan[0] - span[1])/self.dpi)
            else:
                spaces.append((span[0] - pspan[1])/self.dpi)
        # Get the new spacing, remember ratios in like row/column are identical
        # warnings.warn('All panels in this row or column are invisible. You are probably doing something really weird.')
        idx = (-2 if side in 'br' else 1)
        if flush or not spaces:
            space = 0
        else:
            space = max([0, ratios[0][idx] - min(spaces) + pad])
        for iratios in ratios:
            iratios[idx] = space
        # Get array of panel widths, which we left unchanged
        # Then return the 'hpanels' and 'wpanels' arguments
        width = pratios[0][::2]
        return sum(sep) + sum(width) + space

    def smart_tight_layout(self, renderer=None):
        """
        Conforms figure edges to tight bounding box around content without
        messing up subplot aspect ratios and panel widths, and adjusts inner
        spaces so that content from individual subplots (e.g.  axis tick
        labels) does not overlap. This is called automatically when
        `~Figure.draw` or `~Figure.save` are called, unless ``tight=False``
        was passed to `~proplot.subplots.subplots`.

        Parameters
        ----------
        renderer : `~matplotlib.backend_bases.RendererBase`
            The backend renderer. If ``None``, it is retrieved with
            `~matplotlib.figure.Figure.canvas.get_renderer`.
        """
        # Initial stuff
        gridspec = self._main_gridspec
        subplots_kw = self._subplots_kw
        if subplots_kw is None or gridspec is None:
            return

        # Tick fudge factor
        # NOTE: Matplotlib majorTicks, minorTicks attributes contain *un-updated*
        # lists of ticks. To internally update the list and get actual ticks
        # in light of axis locator, use get_major_ticks and get_minor_ticks.
        # WARNING: Could not figure out consistent recipe for calculating tick
        # pad errors in the seemingly infinite number of possible combos of
        # major and minor tick lengths, sides, label sides, and axis label sides.
        # New recommendation is to just use the 'wflush' and 'hflush' overrides
        # if user needs axes directly flush against each other, and in other
        # scenarios error won't matter much -- user can fudge subplotpad if necessary.
        # for ax in self._main_axes:
        #     iaxs = [ax, ax.leftpanel, ax.rightpanel, ax.bottompanel, ax.toppanel]
        #     for iax in iaxs:
        #         if not isinstance(iax, axes.CartesianAxes) or not iax.visible():
        #             continue
        #         # Store error in *points*, because we use it to adjust bounding
        #         # box span, whose default units are in dots.
        #         # Left right ticks
        #         if not iax.yaxis.label.get_text().strip():
        #             side = None
        #         else:
        #             side = iax.yaxis.label_position
        #         # *Always* pad minor tick error, for some reason
        #         pad = (0,0)
        #         ticks = [*iax.yaxis.minorTicks] # copy
        #         ticks_after = iax.yaxis.get_minor_ticks()
        #         if ticks:
        #             lbool, rbool = 1, 1
        #             if ticks_after: # only fix the side where ticks not actually present
        #                 lbool = int(not ticks_after[0].tick1On)
        #                 rbool = int(not ticks_after[0].tick2On)
        #             pad = 1 + ticks[0].get_tick_padding()
        #             pad = (lbool*pad*(side != 'left'), rbool*pad*(side != 'right'))
        #         iax._ytick_pad_error += np.array(pad)*self.dpi/72 # is left, right tuple
        #         # Pad major tick error only if ticks present, and label side
        #         # pad = (0,0)
        #         # ticks = iax.yaxis.majorTicks
        #         # if ticks:
        #         #     pad = ticks[0].get_tick_padding()
        #         #     pad = (pad*ticks[0].tick1On*(side != 'left'),
        #         #            pad*ticks[0].tick2On*(side != 'right'))
        #         # iax._ytick_pad_error += np.array(pad)*self.dpi/72 # is left, right tuple

        #----------------------------------------------------------------------#
        # Put tight box *around* figure
        #----------------------------------------------------------------------#
        if self._smart_tight_outer:
            # Get old un-updated bounding box
            pad = self._smart_borderpad
            obbox = self.bbox_inches # original bbox
            oxmax, oymax = obbox.xmax, obbox.ymax
            # Get new bounding box
            if not renderer: # cannot use the below on figure save! figure becomes a special FigurePDF class or something
                renderer = self.canvas.get_renderer()
            bbox = self.get_tightbbox(renderer)
            xmin, xmax, ymin, ymax = bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax
            # Calculate new edges
            if any(coord is None or np.isnan(coord) for coord in (xmin, xmax, ymin, ymax)):
                warnings.warn('Bounding box has NaNs, cannot get outer tight layout.')
            else:
                loff, boff = xmin, ymin # left bottom margins
                roff, toff = oxmax - xmax, oymax - ymax # top right margin *deltas*
                for key,off in zip(
                    ('left','right','bottom','top'),
                    (loff,roff,boff,toff)
                    ):
                    margin = subplots_kw[key] - off + pad
                    if margin<0:
                        warnings.warn(f'Got negative {key} margin in smart tight layout.')
                    subplots_kw[key] = margin

        #----------------------------------------------------------------------#
        # Prevent overlapping axis tick labels and whatnot *within* figure
        #----------------------------------------------------------------------#
        if self._smart_tight_subplot and self._main_axes:
            # Get bounding box for each axes
            self._tight_bboxs(renderer) # can use same bboxs

            # First for the main axes, then add in the panels
            pad = self._smart_subplotpad
            xspans, yspans, xrange, yrange = _ax_props(self._main_axes, renderer)
            for side,paxs in zip('lbr', (self.leftpanel, self.bottompanel, self.rightpanel)):
                # Initial stuff
                # TODO: Shouldn't this also test if get_visible is on?
                if not any(paxs):
                    continue
                npanel = len(paxs)
                order = paxs._order
                n = paxs._n # rows or columns
                m = npanel//n # columns or rows
                # Add bbox spans and "ranges" for figure panels
                if (side=='b' and order=='C') or (side=='r' and order!='C'):
                    paxs_main = paxs[:n]
                elif (side=='l' and order=='C'):
                    paxs_main = paxs[n-1::n]
                elif (side=='l' and order!='C'):
                    paxs_main = paxs[-n:]
                elif (side=='r' and order=='C') or (side=='b' and order!='C'):
                    paxs_main = paxs[::n]
                xs, ys, xr, yr = _ax_props(paxs_main, renderer)
                xspans, yspans = np.vstack((xspans, xs)), np.vstack((yspans, ys))
                xrange, yrange = np.vstack((xrange, xr)), np.vstack((yrange, yr))
                # Space between stacked figure panels
                # First turn 2D list of panels into list of stacked panels!
                if (side=='b' and order=='C') or (side!='b' and order!='C'):
                    paxs = [paxs[i::n] for i in range(n)] # n is number of panels
                else:
                    paxs = [paxs[i*n:i*n+n] for i in range(m)] # m is number of panels
                sep = self._panel_tight_layout(side, paxs, renderer, figure=True)
                subplots_kw[side + 'sep'] = sep

            # Get *original* "wspace" and "hspace" for comparative purposes
            wspace_orig, hspace_orig = subplots_kw['wspace'], subplots_kw['hspace'] # originals
            if self.leftpanel:
                wspace_orig = [subplots_kw['lspace'], *wspace_orig]
            if self.rightpanel:
                wspace_orig = [*wspace_orig, subplots_kw['rspace']]
            if self.bottompanel:
                hspace_orig = [*hspace_orig, subplots_kw['bspace']]

            # Find groups of axes with touching left/right sides
            # We generate a list (xgroups), each element corresponding to a
            # column of *space*, containing lists of dictionaries with 'l' and
            # 'r' keys, describing groups of "touching" axes in that column
            xgroups = []
            xrange[:,1] += 1
            yrange[:,1] += 1
            ncols = self._main_axes[0]._ncols
            nrows = self._main_axes[0]._nrows
            for cspace in range(1, ncols): # spaces between columns
                groups = []
                for row in range(nrows):
                    filt = (yrange[:,0] <= row) & (row < yrange[:,1])
                    if sum(filt)<=1: # no interface here
                        continue
                    right, = np.where(filt & (xrange[:,0] == cspace))
                    left,  = np.where(filt & (xrange[:,1] == cspace))
                    if not (left.size==1 and right.size==1):
                        continue # normally both zero, but one can be non-zero e.g. if have left and bottom panel, empty space in corner
                    left, right = left[0], right[0]
                    added = False
                    for group in groups:
                        if left in group['l'] or right in group['r']:
                            group['l'].update((left,))
                            group['r'].update((right,))
                            added = True
                            break
                    if not added:
                        groups.append({'l':{left}, 'r':{right}}) # form new group
                xgroups.append(groups)

            # Find groups of axes with touching bottom/top sides
            ygroups = []
            for rspace in range(1, nrows): # spaces between rows
                groups = []
                for col in range(ncols):
                    filt = (xrange[:,0] <= col) & (col < xrange[:,1])
                    top,    = np.where(filt & (yrange[:,0] == rspace)) # gridspec indices go top-to-bottom, left-to-right
                    bottom, = np.where(filt & (yrange[:,1] == rspace))
                    if not (bottom.size==1 and top.size==1):
                        continue
                    bottom, top = bottom[0], top[0]
                    added = False
                    for group in groups:
                        if bottom in group['b'] or top in group['t']:
                            group['b'].update((bottom,))
                            group['t'].update((top,))
                            added = True
                            break
                    if not added:
                        groups.append({'b':{bottom}, 't':{top}}) # form new group
                ygroups.append(groups)

            # Correct wspace and hspace for each group in the 'groups' list
            # For now, ignore the wflush and hflush settings
            wspace = []
            for space_orig, groups in zip(wspace_orig, xgroups):
                if not groups:
                    wspace.append(space_orig)
                    continue
                seps = [] # will use *maximum* necessary separation
                for group in groups: # groups with touching edges
                    left = max(xspans[idx,1] for idx in group['l']) # axes on left side of column
                    right = min(xspans[idx,0] for idx in group['r']) # axes on right side of column
                    seps.append((right - left)/self.dpi)
                wspace.append(max((0, space_orig - min(seps) + pad)))
            hspace = []
            for space_orig, groups in zip(hspace_orig, ygroups):
                if not groups:
                    hspace.append(space_orig)
                    continue
                seps = [] # will use *maximum* necessary separation
                for group in groups: # groups with touching edges
                    bottom = min(yspans[idx,0] for idx in group['b'])
                    top = max(yspans[idx,1] for idx in group['t'])
                    seps.append((bottom - top)/self.dpi)
                hspace.append(max((0, space_orig - min(seps) + pad)))

            # If had figure panels, need to pull out args
            # If wflush or hflush enabled, overwrite the resulting wspace and hspace
            lspace, rspace, bspace = 0, 0, 0 # does not matter if no panels
            if self.leftpanel:
                lspace, *wspace = wspace
            if self.rightpanel:
                *wspace, rspace = wspace
            if self.bottompanel:
                *hspace, bspace = hspace
            if self._subplot_wflush:
                wspace = [0]*len(wspace)
            if self._subplot_hflush:
                hspace = [0]*len(hspace)
            subplots_kw.update({
                'wspace':wspace, 'hspace':hspace,
                'lspace':lspace, 'rspace':rspace, 'bspace':bspace,
                })

        #----------------------------------------------------------------------#
        # The same, but for spaces between *axes panels*
        #----------------------------------------------------------------------#
        # NOTE: any([]) is False so if _main_axes is empty, this is skipped
        panels = False
        if self._smart_tight_panel and any(ax._panels_main_gridspec for ax in self._main_axes):
            # Update bboxs
            panels = True
            for ax in self._main_axes:
                ax.rowlabel.set_visible(False) # needed for panel tight layout!
                ax.collabel.set_visible(False)
            self._tight_bboxs(renderer) # can use same bboxs

            # Bottom, top panels in same rows
            axs = self._main_axes
            ncols = axs[0]._ncols
            nrows = axs[0]._nrows
            hpanels = []
            for row in range(nrows):
                ihpanels = 0
                paxs = [ax.bottompanel for ax in axs if ax._yrange[1]==row]
                ihpanels += self._panel_tight_layout('b', paxs, renderer)
                paxs = [ax.toppanel for ax in axs if ax._yrange[0]==row]
                ihpanels += self._panel_tight_layout('t', paxs, renderer)
                hpanels.append(ihpanels)
            # Left, right panels in same columns
            wpanels = []
            for col in range(ncols):
                iwpanels = 0
                paxs = [ax.leftpanel for ax in axs if ax._xrange[0]==col]
                iwpanels += self._panel_tight_layout('l', paxs, renderer)
                paxs = [ax.rightpanel for ax in axs if ax._xrange[1]==col]
                iwpanels += self._panel_tight_layout('r', paxs, renderer)
                wpanels.append(iwpanels)
            # Adjust for panels
            if self.leftpanel:
                _, *wpanels = wpanels
            if self.rightpanel:
                *wpanels, _ = wpanels
            if self.bottompanel:
                *hpanels, _ = hpanels
            # Add to dictionary
            subplots_kw.update({'wpanels':wpanels, 'hpanels':hpanels})

            # Put back row and column labels
            for ax in self._main_axes:
                ax.rowlabel.set_visible(True) # needed for panel tight layout!
                ax.collabel.set_visible(True)

        #----------------------------------------------------------------------#
        # Update gridspec(s)
        #----------------------------------------------------------------------#
        # Parse arguments and update gridspec
        # Get the new width and height ratios including spaces
        figsize, gridspec_kw, _ = _subplots_kwargs(**subplots_kw)
        gridspec.update(**gridspec_kw)
        self.set_size_inches(figsize)

        if panels:
            # Update width and height ratios of the *inner panel gridspecs* to
            # reflect the new axes heights and widths in the *outer figure gridspec*
            # Note width and height ratios are stored in physical units, inches.
            for ax in self._main_axes:
                # Ratios for outer gridspec
                igridspec = ax._panels_main_gridspec
                if igridspec is None:
                    continue
                idx = (2 if ax.toppanel else 0, 2 if ax.leftpanel else 0)
                wratios = gridspec.get_width_ratios()
                hratios = gridspec.get_height_ratios()
                # Ratios and 'extra' space for axes and its panels
                iwratios = igridspec.get_width_ratios() # gets my *custom* ratios
                ihratios = igridspec.get_height_ratios()
                ihpanels = sum(h for i,h in enumerate(ihratios) if i!=idx[0])
                iwpanels = sum(w for i,w in enumerate(iwratios) if i!=idx[1])
                # Calculate
                xrange = 2*np.array(ax._xrange) # endpoint inclusive
                yrange = 2*np.array(ax._yrange)
                fullwidth = sum(wratios[xrange[0]:xrange[1]+1]) # including panels!
                fullheight = sum(hratios[yrange[0]:yrange[1]+1])
                axwidth = fullwidth - iwpanels
                axheight = fullheight - ihpanels
                # Update axes panels
                iwratios[2 if ax.leftpanel else 0] = axwidth
                ihratios[2 if ax.toppanel else 0] = axheight
                igridspec.set_width_ratios(iwratios)
                igridspec.set_height_ratios(ihratios)
            # Update main gridspec so that it will reflect changes we just made
            # to the panel subplot specs. Gridspec propagate *down*, not *up*.
            gridspec.update(**gridspec_kw)

        # Update axes width and height, used by various functions
        width, height = figsize
        for ax in self._main_axes:
            width_new = abs(ax._position.width)*width
            height_new = abs(ax._position.height)*height
            ax.width, ax.height = width_new, height_new

    def draw(self, renderer, *args, **kwargs):
        """Fixes row and "super" title positions and automatically adjusts the
        main gridspec, then calls the parent `~matplotlib.figure.Figure.draw`
        method."""
        # Prepare for rendering
        self._post_process(renderer)
        # Render
        out = super().draw(renderer, *args, **kwargs)
        return out

    def save(self, *args, **kwargs):
        """Alias for `~Figure.savefig`... because ``Figure.savefig`` is
        completely redundant."""
        return self.savefig(*args, **kwargs)

    def savefig(self, filename, alpha=None, color=None, **kwargs):
        """
        Fixes row and "super" title positions and automatically adjusts the
        main gridspec, then calls the parent `~matplotlib.figure.Figure.savefig`
        method.

        Parameters
        ----------
        filename : str
            The file name and path. Use a tilde ``~`` to represent the home
            directory.
        alpha : float, optional
            Alternative for the `~matplotlib.figure.Figure.savefig`
            `transparent` keyword arg.
        color : color-spec, optional
            Alternative for the `~matplotlib.figure.Figure.savefig`
            `facecolor` keyword arg.
        **kwargs
            Passed to `~matplotlib.figure.Figure.savefig`.
        """
        # Minor changes
        filename = os.path.expanduser(filename)
        if alpha is not None:
            kwargs['transparent'] = not bool(alpha) # 1 is non-transparent
        if color is not None:
            kwargs['facecolor'] = color
            kwargs['transparent'] = False
        # Prepare for rendering
        self._post_process()
        # Render
        return super().savefig(filename, **kwargs) # specify DPI for embedded raster objects

    def add_subplot_and_panels(self, subspec, which=None, order='C', ax_kw={}, *,
            bwidth, bspace, bvisible, bflush, bshare, bsep,
            twidth, tspace, tvisible, tflush, tshare, tsep,
            lwidth, lspace, lvisible, lflush, lshare, lsep,
            rwidth, rspace, rvisible, rflush, rshare, rsep,
            ):
        """
        Creates axes with optional "panels" along the sides of the axes. This
        is called by `subplots` -- you should not have to use this directly.

        Parameters
        ----------
        subspec : `~matplotlib.gridspec.SubplotSpec`
            The `~matplotlib.gridspec.SubplotSpec` instance onto which
            the main subplot and its panels are drawn.
        which : str, optional
            Whether to draw panels on the left, right, bottom, or top
            sides. Should be a string containing any of the characters
            ``'l'``, ``'r'``, ``'b'``, or ``'t'`` in any order.
            Defaults to ``'r'``.
        width, lwidth, rwidth, bwidth, twidth : float or str or list thereof, optional
            Width of left, right, bottom, and top panels, respectively.
            If float or str, widths are same for all panels in the stack. If
            list thereof, specifies widths of each panel in the stack.
            If float, units are inches. If string, units are interpreted by
            `~proplot.utils.units`. Use `width` to set for all sides at once.
        space, lspace, rspace, bspace, tspace : float or str or list thereof, optional
            Empty space between the main subplot and the left, right, bottom,
            and top panels, respectively. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. Use `space` to
            set for all sides at once.
        share, lshare, rshare, bshare, tshare : bool, optional
            Whether to enable axis sharing between *y* axis of main subplot
            and left, right panels, and between *x* axis of main subplot and
            bottom, top panels. Use `share` to set for all sides at once.
        flush, lflush, rflush, bflush, tflush : bool, optional
            Whether *inner* stacked panel should always be flush against the
            subplot, and *stacked* panels flush against each other.
            This overrides the `~Figure.smart_tight_layout` automatic spacing,
            and potentially overrides the `lspace`, `rspace`, `bspace`, `rspace`,
            `lsep`, `rsep`, `bsep`, and `tsep` spacing. Defaults to ``False``.
            Use `flush` to set for all sides at once.
        visible, lvisible, rvisible, bvisible, tvisible : bool, optional
            Used internally. Whether to make the left, right, bottom, or top
            panel invisible. This helps auto-align rows and columns of subplots
            when they don't all have panels on the same sides.
        stack, lstack, rstack, bstack, tstack : int, optional
            The number of optional "stacked" panels on the left, right, bottom,
            and top sides, respectively. Defaults to ``1``. Use `stack` to
            set for all sides at once.
        sep, lsep, rsep, bsep, tsep : float, str, or list thereof, optional
            The separation between stacked panels. If float, units are inches.
            If string, units are interpreted by `~proplot.utils.units`. Ignored
            if the respecitve `stack` keyword arg is None. Use `sep` to set
            for all sides at once.
        """
        # Which panels
        which = _default(which, 'r')
        if re.sub('[lrbt]', '', which): # i.e. other characters are present
            raise ValueError(f'Whichpanels argument can contain characters l (left), r (right), b (bottom), or t (top), instead got "{which}".')
        # Fill space arrays
        wspace = []
        if 'l' in which:
            wspace = [lspace, *wspace]
        if 'r' in which:
            wspace = [*wspace, rspace]
        hspace = []
        if 't' in which:
            hspace = [tspace, *hspace]
        if 'b' in which:
            hspace = [*hspace, bspace]

        # Fix wspace/hspace in inches, using the Bbox from get_postition
        # on the subspec object to determine physical width of axes to be created
        # Fix widths and heights
        bbox = subspec.get_position(self) # valid since axes not drawn yet
        figwidth, figheight = self.get_size_inches()
        boxheight = abs(bbox.height)*figheight
        boxwidth = abs(bbox.width)*figwidth
        height = boxheight - sum(hspace)
        width = boxwidth - sum(wspace)

        # Figure out hratios/wratios
        # Will enforce (main_width + panel_width)/total_width = 1
        wratios = [width]
        if 'l' in which:
            lextra = sum(lwidth) + sum(lsep)
            wratios = [lextra, wratios[0] - lextra, *wratios[1:]]
        if 'r' in which:
            rextra = sum(rwidth) + sum(rsep)
            wratios = [*wratios[:-1], wratios[-1] - rextra, rextra]
        hratios = [height]
        if 't' in which:
            textra = sum(twidth) + sum(tsep)
            hratios = [textra, hratios[0] - textra, *hratios[1:]]
        if 'b' in which:
            bextra = sum(bwidth) + sum(bsep)
            hratios = [*hratios[:-1], hratios[-1] - bextra, bextra]
        if any(w < 0 for w in wratios):
            raise ValueError(f'Left and/or right panel widths too large for available space {width}in.')
        if any(h < 0 for h in hratios):
            raise ValueError(f'Top and/or bottom panel widths too large for available space {height}in.')

        # Make subplotspec
        nrows = 1 + len(re.sub('[^bt]', '', which))
        ncols = 1 + len(re.sub('[^lr]', '', which))
        rows = np.array([w for w in ('t', 'c', 'b') if w in ('c', *which)])
        cols = np.array([w for w in ('l', 'c', 'r') if w in ('c', *which)])
        gs = FlexibleGridSpecFromSubplotSpec(subplot_spec=subspec,
            nrows=nrows, ncols=ncols,
            wspace=wspace, hspace=hspace,
            width_ratios=wratios, height_ratios=hratios)
        # Draw main axes
        r, = np.where('c'==rows)
        c, = np.where('c'==cols)
        ax = self.add_subplot(gs[r[0],c[0]], **ax_kw)
        ax._panels_main_gridspec = gs
        # Draw panels
        for side,width,sep,flush,share,visible, in zip('lrbt',
                (lwidth, rwidth, bwidth, twidth),
                (lsep, rsep, bsep, tsep),
                (lflush, rflush, bflush, tflush),
                (lshare, rshare, bshare, tshare),
                (lvisible, rvisible, bvisible, tvisible),
                ):
            if side not in which:
                continue
            # Settings, and apply flush
            stack = len(width)
            if flush:
                sep = sep*0.0
            if side in 'lr':
                r, = np.where('c'==rows)
                c, = np.where(side==cols)
                wspace, hspace = sep, []
                wratios, hratios = width, [1]
                nrows, ncols = 1, stack
            else:
                r, = np.where(side==rows)
                c, = np.where('c'==cols)
                wspace, hspace = [], sep
                wratios, hratios = [1], width
                nrows, ncols = stack, 1
            # Make gridspec and draw panels
            paxs = []
            name = {'b':'bottom', 't':'top', 'l':'left', 'r':'right'}[side]
            igs = FlexibleGridSpecFromSubplotSpec(
                subplot_spec=gs[r[0],c[0]],
                nrows=nrows, ncols=ncols,
                wspace=wspace, hspace=hspace,
                width_ratios=wratios, height_ratios=hratios)
            for i in range(stack):
                pax = self.add_subplot(igs[i], side=name,
                    share=(visible and share), flush=flush, visible=visible,
                    parent=ax, projection='panel')
                pax._panels_main_gridspec = gs
                pax._panels_stack_gridspec = igs
                paxs += [pax]
            # Add as axes_grid. Support 2D indexing, even though these are
            # always vector stacks, because consistency. See axes_grid docs.
            n = 1 if (side in 'tb' and order=='C') or (side in 'lr' and order!='C') else stack
            paxs = axes_grid(paxs, n=n, order=order)
            setattr(ax, name + 'panel', paxs)
        # Set up axis sharing
        # * Sharex and sharey methods should be called on the 'child' axes,
        #   with labelling preferred on the axes in the first argument.
        # * This block must come *after* the above; sharex_setup and
        #   sharey_setup will detect invisible panels and not share with them.
        # First the x-axis
        bottom = None
        if 'b' in which and bvisible and bshare:
            bottom = ax.bottompanel[-1]
            for iax in (ax, *ax.bottompanel[:-1]):
                iax._sharex_setup(bottom, 3) # parent is *bottom-most* panel
        if 't' in which and tvisible and tshare:
            bottom = bottom or ax
            for iax in ax.toppanel:
                iax._sharex_setup(bottom, 3)
        # Now the y-axis
        left = None
        if 'l' in which and lvisible and lshare:
            left = ax.leftpanel[0]
            for iax in (*ax.leftpanel[1:], ax):
                iax._sharey_setup(left, 3) # parent is *bottom-most* panel
        if 'r' in which and rvisible and rshare:
            left = left or ax
            for iax in ax.rightpanel:
                iax._sharey_setup(left, 3)
        return ax

#-------------------------------------------------------------------------------
# Primary plotting function, used to create figure/axes
#-------------------------------------------------------------------------------
def _panels_sync(side, nums, panels_kw):
    """Makes sure we leave space for invisible panels for subplots *without*
    panels that reside in the same row or column as subplots *with* panels."""
    # Get dict of panels_kw dictionaries in the column or row, and dict
    # of panels_kw dictionaries in the column or row with *active panels*
    all_kw = {num:panels_kw[num] for num in nums}
    on_kw = {num:kw for num,kw in all_kw.items() if side in kw['which']}
    if not on_kw:
        return 0
    # Compare all other panel dicts to the first dict with active panels
    # on the requested side side
    ref_kw = (*on_kw.values(),)[0]
    nsep, nspace, nwidth, nflush = side + 'sep', side + 'space', side + 'width', side + 'flush'
    sep, space, width, flush = ref_kw[nsep], ref_kw[nspace], ref_kw[nwidth], ref_kw[nflush]
    for num,kw in all_kw.items():
        # Make sure stacked-panel separation and panel widths always match
        kw[nsep] = sep
        kw[nwidth] = width
        kw[nflush] = flush
        kw[nspace] = space
        # Match subplot-panel spacing
        if num not in on_kw:
            kw['which'] += side
            kw[side + 'visible'] = False
    # Return the space we have alotted for the panel(s) in this row/column
    return sum(sep) + sum(width) + space

def _panels_kwargs(
    panels, colorbars, legends,
    panels_kw, colorbars_kw=None, legends_kw=None,
    figure=False, ncols=None, nrows=None
    ):
    """Returns standardized keyword args for axes and figure panels."""
    # Get which panels
    kwout = {}
    translate = {'bottom':'b', 'top':'t', 'left':'l', 'right':'r'}
    panels = translate.get(panels, panels)
    legends = translate.get(legends, legends)
    colorbars = translate.get(colorbars, colorbars)
    allpanels = panels + colorbars + legends
    allsides = 'lrb' if figure else 'lrbt'
    if len({*allpanels}) != len(allpanels):
        raise ValueError('You requested the same side for a panel, colorbar, and/or legend.')

    # Fill non-panels with empty args, copy over extra args to potentially
    # raise errors down the line.
    # NOTE: This is done to require keyword-only arguments to _subplots_kwargs
    # and add_subplot_and_panels, so we don't have to duplicate the rc stuff.
    names = ('width', 'sep', 'flush', 'share', 'space')
    if figure:
        names = (*names, 'span')
    else:
        names = (*names, 'visible')
    for side in {*allsides} - {*allpanels}: # add in specific ones
        for name in names:
            kwout[side + name] = None
    extra = (*names, 'stack') # if keyword arg matches this name, ignore it
    regex = re.compile(f'^[tlrb]?({"|".join(extra)})$')
    if colorbars_kw is None:
        colorbars_kw = panels_kw
    if legends_kw is None:
        legends_kw = panels_kw
    for kwargs in (panels_kw, colorbars_kw, legends_kw):
        for key,value in kwargs.items():
            if not regex.match(key):
                kwout[key] = value

    # Define helper function
    def _panel_prop(side, name, defaults):
        """Returns property from the appropriate dictionary, or returns the
        default from the rc.subplots category."""
        if not isinstance(defaults, tuple):
            defaults = 3*(defaults,)
        for check,kwargs,default in zip((panels, colorbars, legends), (panels_kw, colorbars_kw, legends_kw), defaults):
            if side not in check:
                continue
            if isinstance(default, str):
                default = rc['subplots.' + default]
            return _default(kwargs.get(name, None), kwargs.get(side + name, None), default)
    # Get panel widths, account for stacked panels
    # NOTE: Accounts for stacked panels
    for side in allpanels:
        # Simple props
        stack = _panel_prop(side, 'stack', 1)
        share = _panel_prop(side, 'share', (True,False,False))
        kwout[side + 'share'] = share
        # Widths
        width = _panel_prop(side, 'width', ('panelwidth', 'cbarwidth', 'legwidth'))
        width = np.atleast_1d(units(width))
        if len(width)==1:
            width = np.repeat(width, (stack,))
        if len(width)!=stack:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(width)} widths.')
        kwout[side + 'width'] = width
        # Panel separation
        flush = _panel_prop(side, 'flush', False)
        if stack==1 or flush:
            sep = 0
        else:
            default = 'nolabspace' if share else 'ylabspace' if side in 'lr' else 'xlabspace'
            sep = _panel_prop(side, 'sep', (default,default,default))
        sep = np.atleast_1d(units(sep))
        if len(sep)==1:
            sep = np.repeat(sep, (stack-1,))
        if len(sep)!=stack-1:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(sep)} separations.')
        kwout[side + 'sep'] = sep
        kwout[side + 'flush'] = flush

    # Properties for figure panels or axes panels
    if figure:
        for side in allpanels:
            # Space between panels and main subplots grid
            space = _panel_prop(side, 'space', 'xlabspace' if side=='b' else 'ylabspace' if side=='l' else 'nolabspace')
            kwout[side + 'space'] = units(space)
            # Spanning of panels along subplot rows and columns
            nmax = ncols if side=='b' else nrows
            span = _panel_prop(side, 'span', False)
            if np.iterable(span):
                nums = [*span]
                if len(span)!=nmax:
                    raise ValueError(f'Expected {nmax} {side}panel entries, got {len(span)}.')
            elif span:
                nums = [1]*nmax
            else:
                nums = [*range(1,nmax+1)]
            kwout[side + 'span'] = nums
    else:
        # Panel visibility, toggling
        kwout['which'] = allpanels
        for side in allpanels:
            # Space between panels and parent subplot
            if kwout[side + 'share']:
                prop = 'panelspace'
            else:
                prop = ('xlabspace' if side=='b' else 'ylabspace' if side=='l' else 'panelspace')
            space = _panel_prop(side, 'space', prop)
            kwout[side + 'space'] = units(space)
            # Visibility
            kwout[side + 'visible'] = _panel_prop(side, 'visible', True)

    # Return the dictionary
    return kwout

def _subplots_kwargs(nrows, ncols, aspect, xref, yref, *, # ref is the reference axes used for scaling things
    # Args filled with rc settings by main body of subplots()
    wpanels, hpanels,
    width,  height, axwidth, axheight,
    hspace, wspace, hratios, wratios,
    left,   bottom, right,   top,
    # Args filled with rc settings by _panels_kwargs()
    bspan, bwidth, bspace, bsep, bflush, bshare,
    lspan, lwidth, lspace, lsep, lflush, lshare,
    rspan, rwidth, rspace, rsep, rflush, rshare,
    ):
    """Handle complex keyword args and aliases thereof, and determine figure
    sizes such that aspect ratio of first subplot is preserved. Note
    bwidth, bspace, etc. must be supplied or will get error, but this should
    be taken care of by _parse_panels."""
    # Necessary arguments to reconstruct this grid, with defaults filled in
    subplots_kw = {
        'nrows': nrows, 'ncols': ncols, 'aspect': aspect, 'xref': xref, 'yref': yref,
        'width': width, 'height': height, 'axwidth': axwidth, 'axheight': axheight,
        'wpanels': wpanels, 'hpanels': hpanels, 'hspace': hspace, 'wspace': wspace, 'hratios': hratios, 'wratios': wratios,
        'left': left, 'bottom': bottom, 'right': right, 'top': top,
        'bspan': bspan, 'lspan': lspan, 'rspan': rspan,
        'bwidth': bwidth, 'bsep': bsep, 'bspace': bspace, 'bflush': bflush, 'bshare': bshare, # separation between panels
        'rwidth': rwidth, 'rsep': rsep, 'rspace': rspace, 'rflush': rflush, 'rshare': rshare,
        'lwidth': lwidth, 'lsep': lsep, 'lspace': lspace, 'lflush': lflush, 'lshare': lshare,
        }

    # Determine figure size
    # If width and height are not fixed, will scale them to preserve aspect
    # ratio of the first plot
    dx, dy = xref[1]-xref[0], yref[1]-yref[0]
    rwspace = sum(wspace[xref[0]:xref[1]-1])
    rhspace = sum(hspace[yref[0]:yref[1]-1])
    rwratio = ncols*(sum(wratios[slice(*xref)])/dx)/sum(wratios)
    rhratio = nrows*(sum(hratios[slice(*yref)])/dy)/sum(hratios)
    auto_both = (width is None and height is None)
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    auto_neither = (width is not None and height is not None)
    bpanel_space = sum(bwidth) + sum(bsep) + bspace if bspan else 0
    rpanel_space = sum(rwidth) + sum(rsep) + rspace if rspan else 0
    lpanel_space = sum(lwidth) + sum(lsep) + lspace if lspan else 0
    if np.iterable(aspect):
        aspect = aspect[0]/aspect[1]
    # Determine figure dims from axes width and height
    # Default behavior: axes average 2.0 inches wide
    # aspect = aspect*(hratios[0]/np.mean(hratios))
    # NOTE: Account for 'whitespace' that is spanned by our reference axes.
    # When user passes 'axwidth' or 'axheight' they mean the width or height
    # of the entire thing including the part spanning over empty space.
    if auto_both: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = units(rc['subplots.axwidth'])
        if axheight is not None:
            auto_width = True
            axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows # must omit space covered by reference axes
            height = axheight_all + top + bottom + sum(hspace) + sum(hpanels) + bpanel_space
        if axwidth is not None:
            auto_height = True
            axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
            width = axwidth_all + left + right + sum(wspace) + sum(wpanels) + rpanel_space + lpanel_space
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
    # Determine axes widths and heights from requested figure dims
    # Make sure to exclude the 'wpanels' and 'hpanels' props
    else:
        if auto_width or auto_neither:
            axheight_all = height - top - bottom - sum(hspace) - sum(hpanels) - bpanel_space
            axheight = (axheight_all*dy/nrows) + rhspace # reverse engineered from above
        if auto_height or auto_neither:
            axwidth_all = width - left - right - sum(wspace) - sum(wpanels) - rpanel_space - lpanel_space
            axwidth = (axwidth_all*dx/ncols) + rwspace

    # Automatically scale fig dimensions
    # For e.g. common use case [[1,1,2,2],[0,3,3,0]], make sure we still scale
    # the reference axes like a square even though it occupes two columns of gridspec!
    if auto_width:
        axwidth = axheight*aspect
        axwidth_all = ((axwidth - rwspace)/dx/rwratio)*ncols
        width = axwidth_all + left + right + sum(wspace) + sum(wpanels) + rpanel_space + lpanel_space
    elif auto_height:
        axheight = axwidth/aspect
        axheight_all = ((axheight - rhspace)/dy/rhratio)*nrows
        height = axheight_all + top + bottom + sum(hspace) + sum(hpanels) + bpanel_space
    if axwidth_all<0:
        raise ValueError(f"Not enough room for axes (would have width {axwidth_all}). Try using tight=False, increasing figure width, or decreasing 'left', 'right', or 'wspace' spaces.")
    if axheight_all<0:
        raise ValueError(f"Not enough room for axes (would have height {axheight_all}). Try using tight=False, increasing figure height, or decreasing 'top', 'bottom', or 'hspace' spaces.")

    # Make sure the 'ratios' and 'spaces' are in physical units (we cast the
    # former to physical units), easier then to add stuff as below
    wspace = [*wspace]
    hspace = [*hspace]
    wratios = [*(axwidth_all*wratios/sum(wratios) + wpanels)]
    hratios = [*(axheight_all*hratios/sum(hratios) + hpanels)]
    # Now add outer panels, will enforce constant width.
    nrows += int(bool(bspan))
    ncols += int(bool(rspan)) + int(bool(lspan))
    if bspan:
        hratios = hratios + [sum(bwidth) + sum(bsep)]
        hspace  = hspace  + [bspace]
    if lspan:
        wratios = [sum(lwidth) + sum(lsep)] + wratios
        wspace  = [lspace] + wspace
    if rspan:
        wratios = wratios + [sum(rwidth) + sum(rsep)]
        wspace  = wspace  + [rspace]

    # Keyword args for gridspec class
    bottom = bottom/height
    left   = left/width
    top    = 1 - top/height
    right  = 1 - right/width
    gridspec_kw = {
        'nrows': nrows, 'ncols': ncols,
        'left': left, 'bottom': bottom, 'right': right, 'top': top, # so far no panels allowed here
        'wspace': wspace, 'hspace': hspace, 'width_ratios': wratios, 'height_ratios' : hratios,
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
        flush=False, wflush=None, hflush=None,
        left=None, bottom=None, right=None, top=None, # spaces around edge of main plotting area, in inches
        ecols=None, erows=None, # obsolete?
        tight=None, tightborders=None, tightsubplots=None, tightpanels=None,
        borderpad=None, panelpad=None, subplotpad=None,
        span=None,  spanx=1,  spany=1,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        share=None, sharex=3, sharey=3, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        panel=None, panels=None, legend=None, legends=None, colorbar=None, colorbars=None,
        axpanel=None, axlegend=None, axcolorbar=None,
        axpanel_kw=None, axcolorbar_kw=None, axlegend_kw=None,
        axpanels=None, axlegends=None, axcolorbars=None,
        axpanels_kw=None, axcolorbars_kw=None, axlegends_kw=None,
        basemap=False, proj=None, proj_kw=None,
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
        smaller subplots in the bottom row.

        Integers must range from 1 to the number of plots. For example,
        ``[[1, 4]]`` is invalid.

        ``0`` indicates an empty space. For example, ``[[1, 1, 1], [2, 0, 3]]``
        creates one long subplot in the top row with two subplots in the bottom
        row separated by a space.
    ecols, erows : int or list of int, optional
        List, column numbers (starting from 1) that you want *empty*. Generally
        this is used with `ncols` and `nrows`. With `array`, you can
        just use zeros.

    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    width, height : float or str, optional
        The figure width and height. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units` -- for example,
        `width="10cm"` creates a 10cm wide figure.
    journal : str, optional
        String name corresponding to an academic journal standard that is used
        to control the figure width (and height, if specified). Valid names
        are described in a table below.
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
        Passed to `~proplot.gridspec.FlexibleGridSpecBase`. The height
        and width ratios for the subplot grid. Length of `height_ratios`
        must match the number of rows, and length of `width_ratios` must
        match the number of columns.
    hspace, wspace : float or str or list thereof, optional
        If passed, turns off `tightsubplots`.
        These are passed to `~proplot.gridspec.FlexibleGridSpecBase`, denote
        spacing between each column and row of the grid. If float
        or string, expanded into lists of length ``ncols-1`` (for `wspace`) or
        length ``nrows-1`` (for `hspace`). For each element of the list, if float,
        units are inches. If string, units are interpreted by `~proplot.utils.units`.
    top, bottom, left, right : float or str, optional
        If passed, turns off `tightborders`. These are passed to
        `~proplot.gridspec.FlexibleGridSpecBase`, except `right` and
        `top` now refer to the **margin widths** instead of the *x* and *y*
        coordinates for the right and top of the gridspec grid. If float,
        units are inches. If string, units are interpreted by
        `~proplot.utils.units`.

    sharex, sharey, share : {3, 2, 1, 0}, optional
        The "axis sharing level" for the *x* axis, *y* axis, or both
        axes. Options are as follows:

            0. No axis sharing.
            1. Only draw *axis label* on the leftmost column (*y*) or
               bottommost row (*x*) of subplots. Axis tick labels
               still appear on every subplot.
            2. As in 1, but forces the axis limits to be identical. Axis
               tick labels still appear on every subplot.
            3. As in 2, but only show the *axis tick labels* on the
               leftmost column (*y*) or bottommost row (*x*) of subplots.

        This feature can considerably redundancy in your figure.
    spanx, spany, span : bool or {1, 0}, optional
        Whether to use "spanning" axis labels for the *x* axis, *y* axis, or
        both axes. When ``True`` or ``1``, the axis label for the leftmost
        (*y*) or bottommost (*x*) subplot is **centered** on that column or
        row -- i.e. it "spans" that column or row. For example, this means
        for a 3-row figure, you only need 1 y-axis label instead of 3.

        This feature can considerably redundancy in your figure.

        "Spanning" labels also integrate with "shared" axes. For example,
        for a 3-row, 3-column figure, with ``sharey>1`` and ``spany=1``,
        your figure will have **1** *y*-axis label instead of 9.

    proj : str or dict-like, optional
        The map projection name. If string, applies to all subplots. If
        dictionary, values apply to specific subplots, as with `axpanels`.
        If an axes projection is not specified in the dictionary, that axes
        will be Cartesian.

        For example, with ``plot.subplots(ncols=4, proj={1:'mercator', (2,3):'hammer'})``,
        the leftmost subplot is a Mercator projection, the middle 2 are
        Hammer projections, and the rightmost is a normal Cartesian axes.
    proj_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class on instantiation.
        If dictionary of properties, applies globally. If *dictionary of
        dictionaries* of properties, applies to specific subplots, as with `axpanels`.

        For example, with ``plot.subplots(ncols=2, proj='cyl', proj_kw={1:{'lon_0':0}, 2:{'lon_0':180}})``,
        the left subplot is a Plate Carrée projection centered on the prime meridian,
        and the right subplot is centered on the international dateline.
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.
        If boolean, applies to all subplots. If dictionary, values apply to
        specific subplots, as with `axpanels`.

    panel : str, optional
        Specify which sides of the figure should have a "panel".
        This is great for creating global colorbars and legends.
        String should contain any of the characters ``'l'`` (left panel),
        ``'r'`` (right panel), or ``'b'`` (bottom panel). For example,
        ``'br'`` will draw a right and bottom panel.

        Panel axes are stored as ``leftpanel``, ``rightpanel``, and
        ``bottompanel`` attributes on the figure object. They can also
        be accessed by the attribute aliases ``lpanel``, ``rpanel``, and ``bpanel``.
    panels : str, optional
        Similar to `panel`, but default behavior is to assign a panel
        to *every* row or column of subplots. Individual panels can then
        be accessed with e.g. ``fig.leftpanel[0]``, ``fig.leftpanel[1]``.
    colorbar, legend, colorbars, legends : optional
        Identical to `panel` and `panels`, except the *default* panel width is
        more appropriate for a colorbar or legend. The panel can then be
        "filled" with a colorbar or legend with e.g.
        ``fig.bottompanel.colorbar()`` or ``fig.bottompanel.legend()``.
    larray, rarray, barray : bool or list of int, optional
        Defines how figure panels span rows and columns of subplots.
        Each argument is interpreted as follows:

        * If ``True``, this is the default behavior for the `panel`, `colorbar`,
          and `legend` keyword args. Draws *single* panel spanning *all* columns
          or rows of subplots.
        * If ``False``, this is the default behavior for the `panels`, `colorbars`,
          and `legends` keyword args. Draws *separate* panels *for each* column
          or row of subplots.
        * If list of int, this is interpreted like `array`. The integers specify
          panels that span *arbitrary, contiguous* columns or rows
          of subplots.

        For example, ``plot.suplots(ncols=3, bspan=[1, 2, 2])`` draws a panel
        on the bottom of the first column and spanning the bottom of the right
        2 columns, and ``bspan=[0, 2, 2]`` only draws a panel underneath the
        right 2 columns -- as with `array`, the ``0`` indicates an empty space.
    lspan, rspan, bspan : optional
        Aliases for `larray`, `rarray`, and `barray`.
    lspace, rspace, bspace : float, optional
        If passed, turns off `tightsubplots`.
        Space between the edge of the main subplot grid and the left, right,
        and bottom figure panels. If float, units are inches. If string, units
        are interpreted by `~proplot.utils.units`.
    lwidth, rwidth, bwidth : optional
        See `~Figure.add_subplot_and_panels`, usage is identical.
    lstack, rstack, bstack : optional
        See `~Figure.add_subplot_and_panels`, usage is identical.
    lsep, rsep, bsep : optional
        If passed, turns off `tightsubplots`.
        See `~Figure.add_subplot_and_panels`, usage is identical.
    lflush, rflush, bflush : optional
        As in `~Figure.add_subplot_and_panels`, but only controls whether
        *stacked* panels are flush against each other -- i.e. does not make
        figure panels flush against the main subplots.
    lshare, rshare, bshare : optional
        As in `~Figure.add_subplot_and_panels`, but only controls axis
        sharing between *stacked* panels.

    axpanel, axpanels : str or dict-like, optional
        Specifies which axes should have "panels". Both `axpanel` and `axpanels`
        are acceptable, because it is hard to remember otherwise.
        The argument is interpreted as follows:

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
        ``bottompanel`` and ``toppanel`` attributes on axes objects.
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
        Keyword args passed to `~Figure.add_subplot_and_panels` for panels
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


    The current options for the `journal` keyword argument are as follows.
    Feel free to submit a pull request if you'd like to add additional standards.

    ==============================================  ==============================================================  =========================================================================================================================================================
    Key(s)                                          Size description(s)                                             Organization
    ==============================================  ==============================================================  =========================================================================================================================================================
    ``'pnas1'``, ``'pnas2'``, ``'pnas3'``           1-column, 2-column, landscape page                              `Proceedings of the National Academy of Sciences <http://www.pnas.org/page/authors/submission>`__
    ``'ams1'``, ``'ams2'``, ``'ams3'``, ``'ams4'``  1-column, small 2-column, medium 2-column, full 2-column        `American Meteorological Society <https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/>`__
    ``'agu1'``, ``'agu2'``, ``'agu3'``, ``'agu4'``  1-column, 2-column, full height 1-column, full height 2-column  `American Geophysical Union <https://publications.agu.org/author-resource-center/figures-faq/>`__
    ==============================================  ==============================================================  =========================================================================================================================================================
    """
    # TODO: Generalize axes sharing for right y-axes and top x-axes. Enable a secondary
    # axes sharing mode where we *disable ticklabels and labels*, but *do not
    # use the builtin sharex/sharey API*, suitable for complex map projections.
    # For spanning axes labels, right now only detect **x labels on bottom**
    # and **ylabels on top**. Generalize for all subplot edges.
    #--------------------------------------------------------------------------#
    # Array setup
    #--------------------------------------------------------------------------#
    rc._getitem_mode = 0 # ensure still zero; might be non-zero if had error in 'with context' block
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
    array[array==None] = 0 # use zero for placeholder, otherwise have issues
    try:
        array = array.astype(int)
    except (TypeError,ValueError):
        raise ValueError(f'Invalid subplot array {array}. Must be array of integers from 1 to naxs, with 0 representing empty spaces.')
    if ecols:
        ecols = np.atleast_1d(ecols)
        for col in ecols.flat:
            array[:,col-1] = 0
    if erows:
        erows = np.atleast_1d(erows)
        for row in erows.flat:
            array[row-1,:] = 0
    # Enforce rule
    nums = np.unique(array[array!=0])
    naxs = len(nums)
    if {*nums.flat} != {*range(1, naxs+1)}:
        raise ValueError('Invalid subplot array {array}. Numbers must span integers 1 to naxs (i.e. cannot skip over numbers), with 0 representing empty spaces.')
    nrows = array.shape[0]
    ncols = array.shape[1]

    #--------------------------------------------------------------------------#
    # Panels
    #--------------------------------------------------------------------------#
    # Parse outer panel kwargs, consolidate settings
    for ipanel in (panel,legend,colorbar):
        for side in (ipanel or ''):
            if side + 'array' in kwargs:
                kwargs[side + 'span'] = kwargs.pop(side + 'array')
            elif side + 'span' not in kwargs:
                kwargs[side + 'span'] = True # spans all rows or columns
    panels    = _default(panel, panels, '')
    legends   = _default(legend, legends, '')
    colorbars = _default(colorbar, colorbars, '')
    # Get panel props
    kwargs = _panels_kwargs(panels, colorbars, legends, kwargs,
        ncols=ncols, nrows=nrows, figure=True)
    # Warnings
    axpanels    = _default(axpanel, axpanels, '')
    axcolorbars = _default(axcolorbar, axcolorbars, '')
    axlegends   = _default(axlegend, axlegends, '')
    axpanels_kw    = _default(axpanel_kw, axpanels_kw, {})
    axcolorbars_kw = _default(axcolorbar_kw, axcolorbars_kw, {})
    axlegends_kw   = _default(axlegend_kw, axlegends_kw, {})
    if not axpanels and axpanels_kw:
        warnings.warn(f'Ignoring axpanels keyword args: {axpanels_kw}')
    if not axcolorbars and axcolorbars_kw:
        warnings.warn(f'Ignoring axcolorbars keyword args: {axcolorbars_kw}')
    if not axlegends and axlegends_kw:
        warnings.warn(f'Ignoring axlegends keyword args: {axlegends_kw}')
    # Create dictionaries of panel toggles and settings
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    # TODO: Allow separate settings for separate colorbar, legend, etc. panels
    axpanels    = _axes_dict(naxs, axpanels, kw=False, default='')
    axcolorbars = _axes_dict(naxs, axcolorbars, kw=False, default='')
    axlegends   = _axes_dict(naxs, axlegends, kw=False, default='')
    axpanels_kw    = _axes_dict(naxs, axpanels_kw, kw=True)
    axcolorbars_kw = _axes_dict(naxs, axcolorbars_kw, kw=True)
    axlegends_kw   = _axes_dict(naxs, axlegends_kw, kw=True)

    #--------------------------------------------------------------------------#
    # Default behavior
    #--------------------------------------------------------------------------#
    # TODO: Finer level of customization, have tight figure panels mode and
    # tight axes panel mode? Not too important.
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
        if ('space' in key or 'sep' in key) and value is not None}}
    if args:
        tightsubplots = _default(tightsubplots, False)
        if tightsubplots:
            warnings.warn(_args_string('tightsubplots', args))
    # Turn off tightpanels
    args = {key:value for ikw in (axpanels_kw,axcolorbars_kw,axlegends_kw)
        for jkw in ikw.values() for key,value in jkw.items()
        if ('space' in key or 'sep' in key) and value is not None}
    if args:
        tightpanels = _default(tightpanels, False)
        if tightpanels:
            warnings.warn(_args_string('tightpanels', args))

    #--------------------------------------------------------------------------#
    # Shared and spanning axes, panel syncing
    #--------------------------------------------------------------------------#
    # Figure out rows and columns "spanned" by each axes in list, for
    # axis sharing and axis label spanning settings
    sharex = int(_default(share, sharex))
    sharey = int(_default(share, sharey))
    spanx  = _default(span, spanx)
    spany  = _default(span, spany)
    if sharex not in range(4) or sharey not in range(4):
        raise ValueError('Axis sharing options sharex/sharey can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels).')
    # Get some axes properties, where locations are sorted by axes id.
    # NOTE: These ranges are endpoint exclusive, like a slice object!
    axids = [np.where(array==i) for i in np.sort(np.unique(array)) if i>0] # 0 stands for empty
    idxoff = (0, int(bool(kwargs.get('lspan', None))))
    yrange = idxoff[0] + np.array([[y.min(), y.max()+1] for y,_ in axids]) # yrange is shared columns
    xrange = idxoff[1] + np.array([[x.min(), x.max()+1] for _,x in axids])
    xref = xrange[ref-1,:] # range for reference axes
    yref = yrange[ref-1,:]

    # Shared axes: Generate list of base axes-dependent axes pairs
    # That is, find where the minimum-maximum gridspec extent in 'x' for a
    # given axes matches the minimum-maximum gridspec extent for a base axes.
    xgroups_base, xgroups, grouped = [], [], {*()}
    if sharex:
        for i in range(naxs): # axes now have pseudo-numbers from 0 to naxs-1
            matches = (xrange[i,:]==xrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0] # gives ID number of matching_axes, from 0 to naxs-1
            if i not in grouped and matching_axes.size>1:
                # Find all axes that have the same gridspec 'x' extents
                xgroups += [matching_axes]
                # Get bottom-most axis with shared x; should be single number
                xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
            grouped.update(matching_axes) # bookkeeping; record ids that have been grouped already
    ygroups_base, ygroups, grouped = [], [], {*()}
    # Similar for *y* axis
    if sharey:
        for i in range(naxs):
            matches = (yrange[i,:]==yrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups += [matching_axes]
                ygroups_base += [matching_axes[np.argmin(xrange[matching_axes,0])]] # left-most axis with shared y, for matching_axes
            grouped.update(matching_axes) # bookkeeping; record ids that have been grouped already

    # Interpret panel kwargs and apply defaults
    for num in range(1,naxs+1):
        axpanels_kw[num] = _panels_kwargs(
            axpanels[num], axcolorbars[num], axlegends[num],
            axpanels_kw[num], axcolorbars_kw[num], axlegends_kw[num],
            figure=False)
    # Make sure same panels are present on *all rows and columns*
    # Also sompute the extra space alloted for panels
    hpanels = []
    for row in range(nrows):
        ihpanels = 0
        idxs, = np.where(yrange[:,0]==row) # note ranges are *slices*, not endpoint inclusive
        ihpanels += _panels_sync('t', idxs+1, axpanels_kw)
        idxs, = np.where((yrange[:,1]-1)==row)
        ihpanels += _panels_sync('b', idxs+1, axpanels_kw)
        hpanels.append(ihpanels)
    wpanels = []
    for col in range(ncols):
        iwpanels = 0
        idxs, = np.where(xrange[:,0]==col)
        iwpanels += _panels_sync('l', idxs+1, axpanels_kw)
        idxs, = np.where((xrange[:,1]-1)==col)
        iwpanels += _panels_sync('r', idxs+1, axpanels_kw)
        wpanels.append(iwpanels)
    wpanels, hpanels = np.array(wpanels), np.array(hpanels)

    #--------------------------------------------------------------------------#
    # Get basemap.Basemap or cartopy.CRS instances for map, and
    # override aspect ratio.
    #--------------------------------------------------------------------------#
    # NOTE: Cannot have mutable dict as default arg, because it changes the
    # "default" if user calls function more than once! Swap dicts for None.
    basemap = _axes_dict(naxs, basemap, kw=False, default=False)
    proj    = _axes_dict(naxs, _default(proj, 'cartesian'), kw=False, default='cartesian')
    proj_kw = _axes_dict(naxs, _default(proj_kw, {}), kw=True)
    axes_kw = {num:{} for num in range(1, naxs+1)}  # stores add_subplot arguments
    for num,name in proj.items():
        # The default, my CartesianAxes projection
        if name=='cartesian':
            axes_kw[num]['projection'] = 'cartesian'
        # Builtin matplotlib polar axes, just use my overridden version
        elif name=='polar':
            axes_kw[num]['projection'] = 'polar2'
            if num==ref:
                aspect = 1
        # Custom Basemap and Cartopy axes
        elif name:
            package = 'basemap' if basemap[num] else 'cartopy'
            instance, iaspect, kwproj = projs.Proj(name, basemap=basemap[num], **proj_kw[num])
            if num==ref:
                aspect = iaspect
            axes_kw[num].update({'projection':package, 'map_projection':instance})
            axes_kw[num].update(kwproj)
        else:
            raise ValueError('All projection names should be declared. Wut.')

    #--------------------------------------------------------------------------#
    # Figure architecture
    #--------------------------------------------------------------------------#
    # Figure and/or average axes dimensions
    names, values = (), ()
    if journal:
        figsize = _journals(journal) # if user passed width=<string>, will use that journal size
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
    # NOTE: Ratios are scaled to take physical units in _subplots_kwargs, so
    # user can manually provide hspace and wspace in physical units.
    hspace = np.atleast_1d(units(_default(hspace,
        rc['subplots.titlespace'] + rc['subplots.innerspace'] if sharex==3
        else rc['subplots.xlabspace'] if sharex in (1,2) # space for tick labels and title
        else rc['subplots.titlespace'] + rc['subplots.xlabspace'])))
    wspace = np.atleast_1d(units(_default(wspace,
        rc['subplots.innerspace'] if sharey==3
        else rc['subplots.ylabspace'] - rc['subplots.titlespace'] if sharey in (1,2) # space for tick labels only
        else rc['subplots.ylabspace']
        )))
    wratios = np.atleast_1d(_default(width_ratios, wratios, axwidths, 1))
    hratios = np.atleast_1d(_default(height_ratios, hratios, axheights, 1))
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    if flush or wflush:
        wspace = wspace*0.0
    if flush or hflush:
        hspace = hspace*0.0
    if len(wratios)==1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios)==1:
        hratios = np.repeat(hratios, (nrows,))
    left   = units(_default(left,   rc['subplots.ylabspace']))
    bottom = units(_default(bottom, rc['subplots.xlabspace']))
    right  = units(_default(right,  rc['subplots.nolabspace']))
    top    = units(_default(top,    rc['subplots.titlespace']))
    # Parse arguments, fix dimensions in light of desired aspect ratio
    figsize, gridspec_kw, subplots_kw = _subplots_kwargs(
            nrows, ncols, aspect, xref, yref,
            left=left, right=right, bottom=bottom, top=top,
            width=width, height=height, axwidth=axwidth, axheight=axheight,
            wratios=wratios, hratios=hratios, wspace=wspace, hspace=hspace, wpanels=wpanels, hpanels=hpanels,
            **kwargs)
    # Apply settings and add attributes
    gs = FlexibleGridSpec(**gridspec_kw)

    #--------------------------------------------------------------------------#
    # Create blank figure
    #--------------------------------------------------------------------------#
    fig = plt.figure(FigureClass=Figure, tight=tight, figsize=figsize,
        tightborders=tightborders, tightsubplots=tightsubplots, tightpanels=tightpanels,
        borderpad=borderpad, subplotpad=subplotpad, panelpad=panelpad,
        flush=flush, wflush=wflush, hflush=hflush,
        autoformat=autoformat,
        )
    fig._locked = False
    fig._main_gridspec = gs
    fig._subplots_kw = subplots_kw

    #--------------------------------------------------------------------------#
    # Draw on figure
    #--------------------------------------------------------------------------#
    # Create axes
    axs = naxs*[None] # list of axes
    for idx in range(naxs):
        num = idx + 1
        ax_kw = axes_kw[num]
        ax_kw.update({'number':num, 'spanx':spanx, 'spany':spany})
        ax_kw['number'] = num
        if axpanels_kw[num]['which']: # non-empty
            axs[idx] = fig.add_subplot_and_panels(gs[slice(*yrange[idx,:]), slice(*xrange[idx,:])],
                    order=order, ax_kw=ax_kw, **axpanels_kw[num])
        else:
            axs[idx] = fig.add_subplot(gs[slice(*yrange[idx,:]), slice(*xrange[idx,:])],
                    **ax_kw)
    # Create outer panels
    for side in 'blr':
        # Draw panel from gridspec
        panels = subplots_kw[side + 'span']
        if not panels:
            continue
        paxs = []
        name = {'b':'bottom', 'l':'left', 'r':'right'}[side]
        for num in np.unique(panels).flat:
            # Get subspec, and settings
            # TODO: Allow for different width ratios and stuff
            if num==0:
                continue
            off = idxoff[0] if side in 'lr' else idxoff[1]
            idx, = np.where(panels==num)
            idx = slice(off + min(idx), off + max(idx) + 1)
            flush = kwargs[side + 'flush']
            if side=='r':
                subspec = gs[idx,-1]
                wspace, hspace = kwargs['rsep'], []
                wratios, hratios = kwargs['rwidth'], 1
                nrows, ncols = 1, len(wratios)
            elif side=='l':
                subspec = gs[idx,0]
                wspace, hspace = kwargs['lsep'], []
                wratios, hratios = kwargs['lwidth'], 1
                nrows, ncols = 1, len(wratios)
            elif side=='b':
                subspec = gs[-1,idx]
                wspace, hspace = [], kwargs['bsep']
                wratios, hratios = 1, kwargs['bwidth']
                nrows, ncols = len(hratios), 1
            # Make gridspec for containing the "stack" of panels
            ipaxs = []
            igs = FlexibleGridSpecFromSubplotSpec(subplot_spec=subspec,
                    nrows=nrows, ncols=ncols,
                    wspace=wspace, hspace=hspace,
                    width_ratios=wratios, height_ratios=hratios,
                    )
            for i in range(max((nrows,ncols))):
                ipax = fig.add_subplot(igs[i], projection='panel', side=name, flush=flush)
                ipax._panels_main_gridspec = fig._main_gridspec
                ipax._panels_stack_gridspec = igs
                ipaxs += [ipax]
            paxs += [ipaxs]
        # Sort panel axes into row-major or column-major order
        if (side=='b' and order=='C') or (side in 'lr' and order!='C'):
            paxs = [*zip(*paxs)]
        # Store in axes_grid with support for 2D indexing
        n = len(paxs[0])
        paxs = [ax for ipaxs in paxs for ax in ipaxs]
        paxs = axes_grid(paxs, n=n, order=order)
        setattr(fig, name + 'panel', paxs)

    # Set up axis sharing
    # NOTE: You must loop through twice, separately for x and y sharing. Fails
    # otherwise, don't know why.
    if sharex:
        for idx in range(naxs):
            igroup = np.where([idx in g for g in xgroups])[0]
            if igroup.size!=1:
                continue
            sharex_ax = axs[xgroups_base[igroup[0]]]
            if axs[idx] is sharex_ax:
                continue
            axs[idx]._sharex_setup(sharex_ax, sharex)
    if sharey:
        for idx in range(naxs):
            igroup = np.where([idx in g for g in ygroups])[0] # np.where works on lists
            if igroup.size!=1:
                continue
            sharey_ax = axs[ygroups_base[igroup[0]]]
            if axs[idx] is sharey_ax:
                continue
            axs[idx]._sharey_setup(sharey_ax, sharey)
    for ax in axs: # check axes don't belong to multiple groups, should be impossible unless my code is completely wrong...
        for name,groups in zip(('sharex', 'sharey'), (xgroups, ygroups)):
            if sum(ax in group for group in groups)>1:
                raise ValueError(f'Something went wrong; axis {idx:d} belongs to multiple {name} groups.')

    # Return results
    # Store axes in axes_grid with support for 2D indexing
    fig._main_axes = axs
    fig._ref_num = ref
    fig._locked = True
    return fig, axes_grid(axs, n=(ncols if order=='C' else nrows), order=order)

