#!/usr/bin/env python3
"""
The axes subclasses and method wrappers central to this library.
You should focus on the ``format`` and ``smart_update`` methods:

* `BaseAxes.format`
* `BaseAxes.smart_update`
* `XYAxes.smart_update`
* `MapAxes.smart_update`

`BaseAxes.format` calls the ``smart_update`` methods in turn.
This is your one-stop-shop for changing axes settings like
*x* and *y* axis limits, axis labels, tick locations, tick label
format, grid lines, scales, titles, a-b-c labelling, adding
geographic features, and much more.

.. raw:: html

   <h1>Method wrappers</h1>

Wrappers for ProPlot axes methods are documented in the "functions" table.
Method wrapping is done by dynamically overriding the
`~proplot.axes.BaseAxes.__getattribute__` methods on
`~proplot.axes.BaseAxes` and its subclasses.
You may be wondering: why did I choose to do this, rather than using
**decorators**? Two reasons:

1. Brevity. For example: `~proplot.wrappers.cmap_wrapper` overrides **a dozen** different
   methods. This lets me override these methods in *one* line, instead of 50
   lines. To see which methods are overriden, the user can simply check the
   documentation.
2. Documentation. If I wrapped every method, the sphinx `autodoc <http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_
   documentation generator would inherit docstrings from the parent methods.
   In other words, the plotting method docstrings would get duplicated on
   my website from the matplotlib website, generally with a bunch of errors.

   I could also override these methods with my own docstrings, but that would
   mean when the user tries e.g. ``help(ax.contourf)``, they would see my
   own, very brief docstring instead of the comprehensive matplotlib reference
   they were probably looking for! I only add a handful of features to these
   functions, but the native methods generally have way more options.

It should be noted that dynamically wrapping every time the user accesses
the corresponding method will be slower than "decoration", which just wraps them
once on declaration of the axes. But this was not found to significantly affect
performance. And anyway, `Premature Optimization is the Root of All Evil
<http://wiki.c2.com/?PrematureOptimization>`__.
"""
# Note that even if not in IPython notebook, io capture output still works
# functools.wraps preserves __name__ metadata; see comment:
# https://stackoverflow.com/a/739665/4970632
# import os
# import re
import numpy as np
import warnings
from numbers import Number
from matplotlib.cbook import mplDeprecation
import matplotlib.projections as mproj
import matplotlib.axes as maxes
import matplotlib.dates as mdates
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections

# Local modules, projection sand formatters and stuff
from .rcmod import rc, _rc_names_nodots
from . import utils, projs, axistools, wrappers
from .utils import _default, ic, units
from .gridspec import FlexibleGridSpecFromSubplotSpec

# Aliases for panel names
_aliases = {
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'tpanel': 'toppanel',
    'lpanel': 'leftpanel'
    }

# Silly recursive function, returns a...z...aa...zz...aaa...zzz
_abc_string = 'abcdefghijklmnopqrstuvwxyz'
def _abc(i):
    if i < 26:
        return _abc_string[i]
    else:
        return _abc(i - 26) + _abc_string[i % 26]

# Filter warnings, seems to be necessary before drawing stuff for first time,
# otherwise this has no effect (e.g. if you stick it in a function)
warnings.filterwarnings('ignore', category=mplDeprecation)
# Optionally import mapping toolboxes
# Main conda distro says they are incompatible, so make sure not required!
try:
    from cartopy.mpl.geoaxes import GeoAxes
except ModuleNotFoundError:
    GeoAxes = object

#------------------------------------------------------------------------------#
# Generalized custom axes class
# WARNING: Wanted to bulk wrap methods using __new__ on a *metaclass*, since
# this would just wrap method *once* and not every single time user accesses
# object. More elegant, but __new__ *does not receive inherited methods* (that
# comes later down the line), so we can't wrap them. Anyway overriding
# __getattribute__ is fine, and premature optimiztaion is root of all evil!
#------------------------------------------------------------------------------#
def _update_text(obj, **kwargs):
    """Allows updating new text properties introduced by override."""
    # Attempt update
    try:
        obj.update(kwargs)
        return obj
    except Exception:
        pass
    # Destroy original text instance and get its properties
    obj.set_visible(False)
    text = kwargs.pop('text', obj.get_text())
    if 'color' not in kwargs:
        kwargs['color'] = obj.get_color()
    if 'weight' not in kwargs:
        kwargs['weight'] = obj.get_weight()
    if 'fontsize' not in kwargs:
        kwargs['fontsize'] = obj.get_fontsize()
    # Position
    pos = obj.get_position()
    x, y = kwargs.pop('position', (None,None))
    x, y = kwargs.pop('x', x), kwargs.pop('y', y)
    x, y = _default(x, pos[0]), _default(y, pos[1])
    if np.iterable(x):
        x = x[0]
    if np.iterable(y):
        y = y[0]
    # Return new object
    return obj.axes.text(x, y, text, **kwargs)

class BaseAxes(maxes.Axes):
    """Lowest-level axes subclass. Handles titles and axis
    sharing. Adds several new methods and overrides existing ones."""
    name = 'base'
    """The registered projection name."""
    def __init__(self, *args, number=None,
            sharex=None, sharey=None, spanx=None, spany=None,
            sharex_level=0, sharey_level=0,
            **kwargs):
        """
        Parameters
        ----------
        number : None or int
            The subplot number, used for a-b-c labelling (see
            `~BaseAxes.format`).
        sharex_level, sharey_level : {3, 2, 1, 0}, optional
            The "axis sharing level" for the *x* axis, *y* axis, or both
            axes.
        sharex, sharey : None or `BaseAxes`, optional
            Axes to use for *x* and *y* axis sharing. Should correspond
            to the subplot in the bottommost row, leftmost column.
        spanx, spany : None or `BaseAxes`, optional
            Axes to use for the "spanning" *x* and *y* axis labels. Should
            correspond to the subplot in the leftmost column, bottommost row.

        See also
        --------
        `~proplot.subplots.subplots`,
        `XYAxes`, `CartopyAxes`, `BasemapAxes`,
        `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`
        """
        # Properties
        self.number = number # for abc numbering
        self._hatch = None # background hatching
        self._spanx = spanx # boolean toggles, whether we want to span axes labels
        self._spany = spany
        self._yrange = None # geometry, filled later
        self._xrange = None
        self._nrows = None
        self._ncols = None
        # Ugly but necessary
        self._xrotated = False # whether manual rotation was applied
        self._yrotated = False # change default tick label rotation when datetime labels present, if user did not override
        self._abc_inside = False
        self._title_inside = False # toggle this to figure out whether we need to push 'super title' up
        self._gridliner_on = False # whether cartopy gridliners are enabled
        self._is_map = False # needed by wrappers, which can't import this file
        # Children and related properties
        self.bottompanel = EmptyPanel()
        self.toppanel    = EmptyPanel()
        self.leftpanel   = EmptyPanel()
        self.rightpanel  = EmptyPanel()
        self._tight_bbox = None # save these
        self._zoom = None # the 'zoom lines' for inset zoom-in axes
        self._panel_parent = None
        self._inset_parent = None # filled later
        self._inset_children = [] # arbitrary number of insets possible
        self._colorbar_parent = None
        self._colorbar_child = None # the *actual* axes, with content and whatnot; may be needed for tight subplot stuff
        self._auto_colorbar = {} # stores plot handles for filling with a colorbar, as function of location
        self._auto_legend = {} # as above, but for legend
        self._auto_colorbar_kw = {} # keyword args for automatic colorbar() call
        self._auto_legend_kw = {} # as above, but for automatic legend() call
        self._filled = False # turned off when panels filled with colorbar or legend
        self._alty_child = None
        self._altx_child = None
        self._alty_parent = None
        self._altx_parent = None
        self._dualy_scale = None # for scaling units on opposite side of ax, and holding data limits fixed
        self._dualx_scale = None
        self._panels_main_gridspec = None # filled with gridspec used for axes subplot and its panels
        self._panels_stack_gridspec = None # filled with gridspec used for 'stacked' panels

        # Call parent
        super().__init__(*args, **kwargs)

        # Set up axis sharing, save geometry
        width, height = self.figure.get_size_inches()
        self.width  = np.diff(self._position.intervalx)*width # position is in figure units
        self.height = np.diff(self._position.intervaly)*height
        if isinstance(self, maxes.SubplotBase):
            nrows, ncols, subspec = self._topmost_subspec()
            self._yrange = ((subspec.num1 // ncols) // 2, (subspec.num2 // ncols) // 2)
            self._xrange = ((subspec.num1 % ncols) // 2,  (subspec.num2 % ncols) // 2)
            self._nrows = 1 + nrows // 2 # remember, we have rows masquerading as empty spaces!
            self._ncols = 1 + ncols // 2
        # Axis sharing, title stuff, new text attributes
        self._sharex_setup(sharex, sharex_level)
        self._sharey_setup(sharey, sharey_level)
        self._title_transform = self.title.get_transform() # save in case it changes
        self.abc = self.text(0, 0, '') # position tbd
        self.collabel = self.text(0, 0, '', va='bottom', ha='center', transform=self._title_transform)
        self.rowlabel = self.text(0, 0, '', va='center', ha='right', transform=self.transAxes)
        # Apply custom props
        self.format(mode=1)

    @wrappers._expand_methods_list
    def __getattribute__(self, attr, *args):
        """
        Wraps methods when they are requested by the user. See
        `~proplot.wrappers.text_wrapper`, and `~proplot.wrappers.legend_wrapper`
        for more info.

        Disables the redundant methods `wrappers._disabled_methods`.

        Enables the attribute aliases ``bpanel`` for ``bottompanel``,
        ``tpanel`` for ``toppanel``, ``lpanel`` for ``leftpanel``, and
        ``rpanel`` for ``rightpanel``.
        """
        attr = _aliases.get(attr, attr)
        obj = super().__getattribute__(attr, *args)
        # Disabled methods
        for message,attrs in wrappers._disabled_methods.items():
            if attr in attrs:
                raise NotImplementedError(message.format(attr))
        # Non-plotting overrides
        # WARNING: All plotting overrides must come after cmap and cycle
        # wrappers but before check centers and edges overrides, so set
        # them specifically for each plot
        if attr=='text':
            obj = wrappers._text_wrapper(self, obj)
        elif attr=='legend' and not isinstance(self, PanelAxes):
            obj = wrappers._legend_wrapper(self, obj)
        return obj

    def _topmost_subspec(self):
        """Get the top-level SubplotSpec (i.e. the one encompassed by an
        axes and all its panels, if any are present)."""
        subspec = self.get_subplotspec()
        gridspec = subspec.get_gridspec()
        while isinstance(gridspec, mgridspec.GridSpecFromSubplotSpec):
            try:
                subspec = gridspec._subplot_spec
            except AttributeError:
                raise ValueError('The _subplot_spec attribute is missing from this GridSpecFromSubplotSpec. Cannot determine the parent GridSpec rows/columns occupied by this slot.')
            gridspec = subspec.get_gridspec()
        nrows, ncols = gridspec.get_geometry()
        return nrows, ncols, subspec

    def _share_short_axis(self, share, side, level):
        """When sharing main subplots, shares the short axes of their side panels."""
        if isinstance(self, PanelAxes):
            return
        axis = 'x' if side[0] in 'lr' else 'y'
        paxs1 = getattr(self, side + 'panel') # calling this means, share properties on this axes with input 'share' axes
        paxs2 = getattr(share, side + 'panel')
        if not all(pax and pax.get_visible() and not pax._filled for pax in paxs1) or \
           not all(pax and pax.get_visible() and not pax._filled for pax in paxs2):
            return
        if len(paxs1) != len(paxs2):
            raise RuntimeError('Sync error. Different number of stacked panels along axes on like column/row of figure.')
        for pax1,pax2 in zip(paxs1,paxs2):
            getattr(pax1, '_share' + axis + '_setup')(pax2, level)

    def _share_long_axis(self, share, side, level):
        """When sharing main subplots, shares the long axes of their side panels,
        assuming long axis sharing is enabled for that panel."""
        if isinstance(self, PanelAxes):
            return
        axis = 'x' if side[0] in 'tb' else 'y'
        paxs = getattr(self, side + 'panel') # calling this means, share properties on this axes with input 'share' axes
        if not all(pax and pax.get_visible() and not pax._filled for pax in paxs) or \
           not all(pax._share for pax in paxs):
            return
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share, level)

    def _sharex_setup(self, sharex, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharex is None or sharex is self:
            return
        if isinstance(self, MapAxes) or isinstance(sharex, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels 
        self._share_short_axis(sharex, 'left',   level)
        self._share_short_axis(sharex, 'right',  level)
        self._share_long_axis(sharex,  'bottom', level)
        self._share_long_axis(sharex,  'top',    level)
        # Builtin features
        self._sharex = sharex
        if level>1:
            self._shared_x_axes.join(self, sharex)
        # "Shared" axis and tick labels
        # WARNING: Assigning *another* axis label to this axes will raise error,
        # because matplotlib tries to draw same Artist twice. Just make it invisible.
        if level>2:
            for t in self.xaxis.get_ticklabels():
                t.set_visible(False)
        self.xaxis.label.set_visible(False)

    def _sharey_setup(self, sharey, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharey is None or sharey is self:
            return
        if isinstance(self, MapAxes) or isinstance(sharey, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels
        self._share_short_axis(sharey, 'bottom', level)
        self._share_short_axis(sharey, 'top',    level)
        self._share_long_axis(sharey,  'left',   level)
        self._share_long_axis(sharey,  'right',  level)
        # Builtin features
        self._sharey = sharey
        if level>1:
            self._shared_y_axes.join(self, sharey)
        # "Shared" axis and tick labels
        if level>2:
            for t in self.yaxis.get_ticklabels():
                t.set_visible(False)
        self.yaxis.label.set_visible(False)

    def _parse_title_args(self, abc=False, pos=None, border=None, linewidth=None, **kwargs):
        """Position title text to the left, center, or right and either
        inside or outside the axes (default is center, outside)."""
        # Maybe all this crap was done already in first call to format
        if pos is None:
            return kwargs
        # Determine position
        ypad = rc.get('axes.titlepad')/(72*self.height) # to inches --> to axes relative
        xpad = rc.get('axes.titlepad')/(72*self.width)  # why not use the same for x?
        xpad_i, ypad_i = xpad*1.5, ypad*1.5 # inside labels need a bit more room
        if not isinstance(pos, str):
            # Coordinate position
            ha = va = 'center'
            x, y = pos
            transform = self.transAxes
        else:
            # Get horizontal position
            if not any(c in pos for c in 'lcr'):
                pos += 'c'
            if not any(c in pos for c in 'oi'):
                pos += 'o'
            if 'c' in pos:
                x = 0.5
                ha = 'center'
            elif 'l' in pos:
                x = 0 + xpad_i*('i' in pos)
                ha = 'left'
            elif 'r' in pos:
                x = 1 - xpad_i*('i' in pos)
                ha = 'right'
            # Get vertical position
            if 'o' in pos:
                y = 1 # leave it alone, may be adjusted during draw-time to account for axis label (fails to adjust for tick labels; see notebook)
                va = 'bottom' # matches title alignment maybe
                transform = self._title_transform
                kwargs['border'] = False
            elif 'i' in pos:
                y = 1 - ypad_i
                va = 'top'
                transform = self.transAxes
                if border is not None:
                    kwargs['border'] = border
                if border and linewidth:
                    kwargs['linewidth'] = linewidth
        # Record _title_inside so we can automatically deflect suptitle
        # If *any* object is outside (title or abc), want to deflect it up
        inside = (not isinstance(pos, str) or 'i' in pos)
        if abc:
            self._abc_inside = inside
        else:
            self._title_inside = inside
        return {'x':x, 'y':y, 'ha':ha, 'va':va, 'transform':transform, **kwargs}

    def format(self, *, mode=2, rc_kw=None, **kwargs):
        """
        Sets up temporary rc settings and calls `~BaseAxes.smart_update`.

        Parameters
        ----------
        rc_kw : None or dict, optional
            A dictionary containing "rc" configuration settings that will
            be applied to this axes. Temporarily updates the
            `~proplot.rcmod.rc` object. See `~proplot.rcmod` for details.
        **kwargs
            Any of three options:

            * A keyword arg for `BaseAxes.smart_update`,
              `XYAxes.smart_update`, or `MapAxes.smart_update`.
            * A global "rc" keyword arg, like ``linewidth`` or ``color``.
            * A standard "rc" keyword arg **with the dots omitted**.
              For example, ``land.color`` becomes ``landcolor``.

            The latter two options update the `~proplot.rcmod.rc`
            object, just like `rc_kw`.

        Other parameters
        ----------------
        mode : int, optional
            The "getitem mode". This is used under-the-hood; you shouldn't
            have to use it directly. Determines whether queries to the
            `~proplot.rcmod.rc` object will ignore :ref:`rcParams`.
            This can help prevent a massive number of unnecessary lookups
            when the settings haven't been changed by the user.
            See `~proplot.rcmod.rc_configurator` for details.

        See also
        --------
        `~proplot.rcmod`,
        `BaseAxes.smart_update`, `XYAxes.smart_update`, `MapAxes.smart_update`
        """
        # Figure out which kwargs are valid rc settings
        # WARNING: First part will fail horribly if mode is not zero!
        # WARNING: If declare {} as kwarg, can get filled by user, so default
        # changes! Need to make default None.
        rc_kw = rc_kw or {}
        kw_extra = {}
        for key,value in kwargs.items():
            key_fixed = _rc_names_nodots.get(key, None)
            if key_fixed is None:
                kw_extra[key] = value
            else:
                rc_kw[key_fixed] = value
        rc._getitem_mode = 0 # might still be non-zero if had error
        # Apply custom settings, potentially
        # TODO: Put this somewhere else?
        kw = {}
        if isinstance(self, PanelAxes) and self._flush:
            if self._side=='right':
                kw.update({'yloc':'right', 'ylabelloc':'right', 'yticklabelloc':'right'})
            elif self._side=='top':
                kw.update({'xloc':'top', 'xlabelloc':'top', 'xticklabelloc':'top'})
        kw.update(kw_extra)
        with rc.context(rc_kw, mode=mode):
            self.smart_update(**kw)

    def smart_update(self, title=None, abc=None,
        figtitle=None, suptitle=None, collabels=None, rowlabels=None, # label rows and columns
        top=True, # nopanel optionally puts title and abc label in main axes
        ):
        """
        Function for formatting axes titles and labels.

        Parameters
        ----------
        title : None or str, optional
            The axes title.
        abc : None or bool, optional
            Whether to apply "a-b-c" subplot labelling based on the
            `number` attribute.

            If `number` is >26, labels will loop around to a, ..., z, aa,
            ..., zz, aaa, ..., zzz, ... God help you if you ever need that
            many labels.
        top : bool, optional
            Whether to try to put title and a-b-c label above the top
            axes panel, if it exists (``True``), or to always put them on
            the main subplot (``False``). Defaults to ``True``.
        rowlabels, colllabels : None or list of str, optional
            The subplot row and column labels. If list, length must match
            the number of subplot rows, columns.
        figtitle, suptitle : None or str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.

            This is more sophisticated than matplotlib's builtin "super title",
            which is just centered between the figure edges and offset from
            the top edge.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.rcmod`,
        `XYAxes.smart_update`, `MapAxes.smart_update`
        """
        # Figure patch (for some reason needs to be re-asserted even if
        # declared before figure is drawn)
        # Look into `~matplotlib.axes.SubplotBase.is_last_row` and
        # `~matplotlib.axes.SubplotBase.is_first_column` methods.
        kw = rc.fill({'facecolor':'figure.facecolor'})
        self.figure.patch.update(kw)

        # Super title and labels
        # NOTE: These are actually *figure-wide* settings, but that line seems
        # to get blurred -- where we have shared axes, spanning labels, and
        # whatnot. May result in redundant assignments if formatting more than
        # one axes, but operations are fast so some redundancy is nbd.
        fig = self.figure # the figure
        suptitle = figtitle or suptitle
        if suptitle is not None:
            kw = rc.fill({
                'fontsize': 'suptitle.fontsize',
                'weight':   'suptitle.weight',
                'color':    'suptitle.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._suptitle_setup(suptitle, **kw)
        if rowlabels is not None:
            kw = rc.fill({
                'fontsize': 'rowlabel.fontsize',
                'weight':   'rowlabel.weight',
                'color':    'rowlabel.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._rowlabels(rowlabels, **kw)
        if collabels is not None:
            kw = rc.fill({
                'fontsize': 'collabel.fontsize',
                'weight':   'collabel.weight',
                'color':    'collabel.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._collabels(collabels, **kw)

        # Axes for title or abc
        pax = self.toppanel[0]
        if top and pax and pax.get_visible() and not pax._filled:
            tax = self.toppanel[0]
        else:
            tax = self
        tax = tax._altx_child or tax # always on top!

        # Create axes title
        # TODO: We check _filled here, but no support for filled top
        # panels; maybe add support?
        # NOTE: Aligning title flush against left or right of axes is
        # actually already a matplotlib feature! Use set_title(loc=loc), and
        # it assigns text to a different hidden object. My version is just
        # more flexible, allows specifying arbitrary postiion.
        kw = rc.fill({
            'pos':       'title.pos',
            'border':    'title.border',
            'linewidth': 'title.linewidth',
            'fontsize':  'title.fontsize',
            'weight':    'title.weight',
            'fontname':  'fontname'
            }, cache=True)
        kw = tax._parse_title_args(**kw)
        if title:
            kw['text'] = title
        if kw:
            tax.title = _update_text(tax.title, **kw)

        # Initial text setup
        # Will only occur if user requests change, or on initial run
        abcformat = rc['abc.format']
        if abcformat and self.number is not None:
            if 'a' not in abcformat and 'A' not in abcformat:
                raise ValueError(f'Invalid abcformat "{abcformat}". Must include letter "a" or "A".')
            abcedges = abcformat.split('a' if 'a' in abcformat else 'A')
            text = abcedges[0] + _abc(self.number-1) + abcedges[-1]
            if 'A' in abcformat:
                text = text.upper()
            tax.abc.set_text(text)
        # Apply any changed or new settings
        kw = rc.fill({
            'pos':       'abc.pos',
            'border':    'abc.border',
            'linewidth': 'abc.linewidth',
            'fontsize':  'abc.fontsize',
            'weight':    'abc.weight',
            'color':     'abc.color',
            'fontname':  'fontname'
            }, cache=True)
        kw = tax._parse_title_args(abc=True, **kw)
        if kw:
            tax.abc = _update_text(tax.abc, **kw)
        tax.abc.set_visible(bool(abc))

    def colorbar(self, *args, loc=None, xspace=None, pad=None, width=None,
        length=None, label=None, extendsize=None,
        **kwargs):
        """
        Adds an *inset* colorbar, sort of like `~matplotlib.axes.Axes.legend`.

        Parameters
        ----------
        loc : None or str or int, optional
            The colorbar location. Just like ``loc`` for the native matplotlib
            `~matplotlib.axes.Axes.legend`, but filtered to only corner
            positions. Options are:

            * ``'ur'`` or ``'upper right'``
            * ``'ul'`` or ``'upper left'``
            * ``'ll'`` or ``'lower left'``
            * ``'lr'`` or ``'lower right'``

            Default is from the rc configuration.
        xspace : None or str or float, optional
            Space allocated for the bottom x-label of the colorbar.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        pad : None or str or float, optional
            Space between the axes edge and the colorbar.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        width : None or str or float, optional
            Colorbar width.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        length : None or str or float, optional
            Colorbar length.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        **kwargs, label, extendsize
            Passed to `~proplot.wrappers.colorbar_wrapper`.
        """
        # Default props
        loc = _default(loc, rc.get('colorbar.loc'))
        extend = units(_default(extendsize, rc.get('colorbar.extendinset')))
        length = units(_default(length, rc.get('colorbar.length')))/self.width
        width = units(_default(width, rc.get('colorbar.width')))/self.height
        pad = units(_default(pad, rc.get('colorbar.axespad')))
        xpad = pad/self.width
        ypad = pad/self.height
        # Space for labels
        xspace = units(_default(xspace, rc.get('colorbar.xspace'))) # for x label
        if not label:
            xspace -= 1.2*rc.get('font.size')/72
        xspace /= self.height
        # Get location in axes-relative coordinates
        if loc in ('upper right','ur'):
            bounds = (1-xpad-length, 1-ypad-width, length, width)
        elif loc in ('upper left','ul'):
            bounds = (xpad, 1-ypad-width, length, width) # x0, y0, width, height in axes-relative coordinate to start
        elif loc in ('lower left','ll'):
            bounds = (xpad, ypad+xspace, length, width)
        elif loc in ('lower right','lr','b','best'):
            bounds = (1-xpad-length, ypad+xspace, length, width)
        else:
            raise ValueError(f'Invalid location {loc}.')
        # Test
        bbox = mtransforms.Bbox.from_bounds(*bounds)
        bb = mtransforms.TransformedBbox(bbox, self.transAxes)
        tr = self.figure.transFigure.inverted()
        bb = mtransforms.TransformedBbox(bb, tr)
        # Axes
        locator = self._make_inset_locator(bounds, self.transAxes)
        bb = locator(None, None)
        ax = maxes.Axes(self.figure, bb.bounds, zorder=5)
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        kwargs.update({'ticklocation':'bottom', 'extendsize':extend, 'label':label})
        return wrappers.colorbar_wrapper(ax, *args, **kwargs)

    def cmapline(self, *args, values=None,
            cmap=None, norm=None,
            interp=0, **kwargs):
        """
        Creates a "colormap line", i.e. a line whose color changes as a function
        of the coordinate `values`. This is actually a collection of lines,
        added as a `~matplotlib.collections.LineCollection` instance. See `this
        matplotlib example <https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html>`__.
        This method is invoked if you call `~matplotlib.axes.Axes.plot` with
        the ``cmap`` keyword arg.

        Parameters
        ----------
        values : list of float
            Values corresponding to points on the line.
        cmap : None or colormap spec, optional
            Colormap specifier, passed to `~proplot.colortools.Colormap`.
        norm : None or `~matplotlib.colors.Normalize`, optional
            The normalizer, used for mapping `values` to colormap colors.
        interp : int, optional
            Number of values between each line joint and each *halfway* point
            between line joints to which you want to interpolate.

        """
        # WARNING: So far this only works for 1D *x* and *y* coordinates. Cannot
        # draw multiple colormap lines at once, unlike `~matplotlib.axes.Axes.plot`.
        # First error check
        if values is None:
            raise ValueError('Requires a "values" keyword arg.')
        if len(args) not in (1,2):
            raise ValueError(f'Requires 1-2 arguments, got {len(args)}.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(y.shape[-1]) if len(args)==1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim!=1 or y.ndim!=1 or values.ndim!=1:
            raise ValueError(f'x ({x.ndim}-d), y ({y.ndim}-d), and values ({values.ndim}-d) must be 1-dimensional.')
        if len(x)!=len(y) or len(x)!=len(values) or len(y)!=len(values):
            raise ValueError(f'{len(x)} xs, {len(y)} ys, but {len(values)} colormap values.')

        # Next draw the line
        # Interpolate values to optionally allow for smooth gradations between
        # values (bins=False) or color switchover halfway between points (bins=True)
        # Next optionally interpolate the corresponding colormap values
        # NOTE: We linearly interpolate here, but user might use a normalizer that
        # e.g. performs log before selecting linear color range; don't need to
        # implement that here
        if interp>0:
            xorig, yorig, vorig = x, y, values
            x, y, values = [], [], []
            for j in range(xorig.shape[0]-1):
                idx = (slice(None, -1) if j+1<xorig.shape[0]-1 else slice(None))
                x.extend(np.linspace(xorig[j], xorig[j+1], interp + 2)[idx].flat)
                y.extend(np.linspace(yorig[j], yorig[j+1], interp + 2)[idx].flat)
                values.extend(np.linspace(vorig[j], vorig[j+1], interp + 2)[idx].flat)
            x, y, values = np.array(x), np.array(y), np.array(values)
        coords = []
        levels = utils.edges(values)
        for j in range(y.shape[0]):
            # Get x/y coordinates and values for points to the 'left' and
            # 'right' of each joint. Also prevent duplicates.
            if j==0:
                xleft, yleft = [], []
            else:
                xleft = [(x[j-1] + x[j])/2, x[j]]
                yleft = [(y[j-1] + y[j])/2, y[j]]
            if j+1==y.shape[0]:
                xright, yright = [], []
            else:
                xleft  = xleft[:-1] # prevent repetition when joined with xright/yright
                yleft  = yleft[:-1] # actually need numbers of x/y coordinates to be same for each segment
                xright = [x[j], (x[j+1] + x[j])/2]
                yright = [y[j], (y[j+1] + y[j])/2]
            pleft  = np.stack((xleft,  yleft), axis=1)
            pright = np.stack((xright, yright), axis=1)
            coords.append(np.concatenate((pleft, pright), axis=0))

        # Create LineCollection and update with values
        hs = mcollections.LineCollection(np.array(coords), cmap=cmap, norm=norm,
                linestyles='-', capstyle='butt', joinstyle='miter')
        hs.set_array(np.array(values))
        hs.update({key:value for key,value in kwargs.items() if key not in ('color',)})

        # Add collection, with some custom attributes
        self.add_collection(hs)
        self.autoscale_view() # data limits not updated otherwise
        hs.values = values
        hs.levels = levels # needed for other functions some
        return hs

#------------------------------------------------------------------------------#
# Specific classes, which subclass the base one
#------------------------------------------------------------------------------#
def _rcloc_to_stringloc(xy, string): # figures out string location
    """Gets location string from the boolean rc settings, for a given string
    prefix like ``'axes.spines'`` or ``'xtick'``."""
    if xy=='x':
        top = rc[f'{string}.top']
        bottom = rc[f'{string}.bottom']
        if top is None and bottom is None:
            return None
        elif top and bottom:
            return 'both'
        elif top:
            return 'top'
        elif bottom:
            return 'bottom'
        else:
            return 'neither'
    elif xy=='y':
        left = rc[f'{string}.left']
        right = rc[f'{string}.right']
        if left is None and right is None:
            return None
        elif left and right:
            return 'both'
        elif left:
            return 'left'
        elif right:
            return 'right'
        else:
            return 'neither'
    else:
        raise ValueError(f'"xy" must equal "x" or "y".')

class XYAxes(BaseAxes):
    """
    Axes subclass for ordinary Cartesian axes. Adds several new methods and
    overrides existing ones.

    See also
    --------
    `~proplot.subplots.subplots`, `BaseAxes`, `PanelAxes`
    """
    name = 'xy'
    """The registered projection name."""
    def __init__(self, *args, **kwargs):
        # Create simple x by y subplot.
        super().__init__(*args, **kwargs)
        # Change the default formatter
        # Mine trims trailing zeros, but numbers no longer aligned. Matter
        # of taste really; will see if others like it.
        formatter = axistools.Formatter('default')
        self.xaxis.set_major_formatter(formatter)
        formatter = axistools.Formatter('default')
        self.yaxis.set_major_formatter(formatter)
        # Reset this; otherwise matplotlib won't automatically change
        # formatter when it encounters certain types of data, like datetime.
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True
        # For tight layout; matplotlib draws tick bbox around ticks, even
        # if they are actually all turned off.
        self._xtick_pad_error = (0,0)
        self._ytick_pad_error = (0,0)

    def __getattribute__(self, attr, *args):
        """Wraps methods when they are requested by the user. See
        `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`,
        `~proplot.wrappers.check_centers`, and `~proplot.wrappers.check_edges`
        for more info."""
        obj = super().__getattribute__(attr, *args)
        # Color wrappers
        if attr in wrappers._cmap_methods: # must come first!
            obj = wrappers._cmap_wrapper(self, obj)
        elif attr in wrappers._cycle_methods:
            obj = wrappers._cycle_wrapper(self, obj)
        # Plotting wrappers
        if attr=='plot':
            obj = wrappers._plot_wrapper(self, obj)
        elif attr=='scatter':
            obj = wrappers._scatter_wrapper(self, obj)
        elif attr in wrappers._centers_methods:
            obj = wrappers._check_centers(self, obj)
        elif attr in wrappers._edges_methods:
            obj = wrappers._check_edges(self, obj)
        return obj

    def smart_update(self,
        xloc=None,          yloc=None,          # aliases for 'where to put spine'
        xmargin=None,       ymargin=None,
        xcolor=None,        ycolor=None,        # separate color for x or y axis spines, ticks, tick labels, and axis labels; useful for twin axes
        xspineloc=None,     yspineloc=None,     # deals with spine options
        xtickloc=None,      ytickloc=None,      # which spines to draw ticks on
        fixticks=False, # whether to always transform locator to FixedLocator
        xlabelloc=None,     ylabelloc=None,
        xticklabelloc=None, yticklabelloc=None, # where to put tick labels
        xtickdir=None,      ytickdir=None,      # which direction ('in', 'out', or 'inout')
        xgrid=None,         ygrid=None,         # gridline toggle
        xgridminor=None,    ygridminor=None,    # minor grids on/off (if ticks off, grid will always be off)
        xtickminor=True,    ytickminor=True,    # minor ticks on/off
        xticklabeldir=None, yticklabeldir=None, # which direction to draw labels
        xtickrange=None,    ytickrange=None,    # limit regions where we assign ticklabels to major-ticks
        xreverse=False,     yreverse=False,     # special properties
        xlabel=None,        ylabel=None,        # axis labels
        xlim=None,          ylim=None,
        xbounds=None,       ybounds=None,       # limit spine bounds?
        xscale=None,        yscale=None,
        xrotation=None,     yrotation=None,     # tick label rotation
        xformatter=None, yformatter=None, xticklabels=None, yticklabels=None,
        xticks=None, xminorticks=None, xlocator=None, xminorlocator=None,
        yticks=None, yminorticks=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
        xlabel_kw={}, ylabel_kw={},
        xscale_kw={}, yscale_kw={},
        xlocator_kw={}, ylocator_kw={},
        xformatter_kw={}, yformatter_kw={},
        xminorlocator_kw={}, yminorlocator_kw={},
        **kwargs): # formatter
        """
        Format the *x* and *y* axis labels, tick locations, tick label
        format, scales, and more.

        Parameters
        ----------
        xlabel, ylabel : None or str, optional
            The *x* and *y* axis labels. Applied with
            `~matplotlib.axes.Axes.set_xlabel`
            and `~matplotlib.axes.Axes.set_ylabel`.
        xlabel_kw, ylabel_kw : dict-like, optional
            The *x* and *y* axis label settings. Applied with the
            `~matplotlib.artist.Artist.update` method on the
            `~matplotlib.text.Text` instance. Options include ``'color'``,
            ``'size'``, and ``'weight'``.
        xlim, ylim : None or length-2 list of float or None, optional
            The *x* and *y* axis data limits. Applied with
            `~matplotlib.axes.Axes.set_xlim` and `~matplotlib.axes.Axes.set_ylim`.
        xreverse, yreverse : bool, optional
            Simply reverses the order of the `xlim` or `ylim` tuples if
            they were provided.
        xbounds, ybounds : None or length-2 list of float, optional
            The *x* and *y* axis data bounds within which to draw the spines.
            For example, this can be used to separate the *x*
            and *y* axis spines, so they don't meet at the corner.
        xtickrange, ytickrange : None or length-2 list of float, optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``[-1,1]`` with
            axis range ``[-5,5]`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1.
        fixticks : bool, optional
            Whether to always transform the tick locators to a
            `~matplotlib.ticker.FixedLocator` instance. Defaults to ``False``.
            If your axis ticks are doing weird things (for example, ticks
            drawn outside of the axis spine), try setting this to ``True``.
        xscale, yscale : None or scale spec, optional
            The *x* and *y* axis scales. Arguments are passed to and interpreted
            by `~proplot.axistools.Scale`. Examples include ``'linear'``
            and ``('cutoff', 0.5, 2)``.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.axistools.Scale`.
        xloc, yloc
            Aliases for `xspineloc`, `yspineloc`.
        xspineloc, yspineloc : {'both', 'bottom', 'top', 'left', 'right', 'neither', 'center', 'zero'}, optional
            The *x* and *y* axis spine locations.
        xtickloc, ytickloc : {'both', 'bottom', 'top', 'left', 'right', 'neither'}, optional
            Which *x* and *y* axis spines should have major and minor tick marks.
        xtickminor, ytickminor : None or bool, optional
            Whether to draw minor ticks on the *x* and *y* axes.
        xtickdir, ytickdir : {'out', 'in', 'inout'}
            Direction that major and minor tick marks point for the *x* and
            *y* axis.
        xgrid, ygrid : None or bool, optional
            Whether to draw major gridlines on the *x* and *y* axis.
        xgridminor, ygridminor : None or bool, optional
            Whether to draw minor gridlines for the *x* and *y* axis.
        xticklabeldir, yticklabeldir : {'in', 'out'}
            Whether to place *x* and *y* axis tick label text inside
            or outside the axes.
        xticks, yticks
            Aliases for `xlocator`, `ylocator`.
        xlocator, ylocator : None or locator spec, optional
            The *x* and *y* axis tick mark positions. Passed to the
            `~proplot.axistools.Locator` constructor.
        xlocator_kw, ylocator_kw : dict-like, optional
            The *x* and *y* axis locator settings. Passed to
            `~proplot.axistools.Locator`.
        xminorticks, yminorticks
            Aliases for `xminorlocator`, `yminorlocator`.
        xminorlocator, yminorlocator
            As for `xlocator`, `ylocator`, but for the minor ticks.
        xminorlocator_kw, yminorlocator_kw
            As for `xlocator_kw`, `ylocator_kw`, but for the minor locator.
        xticklabels, yticklabels
            Aliases for `xformatter`, `yformatter`.
        xformatter, yformatter : None or format spec, optional
            The *x* and *y* axis tick label format. Passed to the
            `~proplot.axistools.Formatter` constructor.
        xformatter_kw, yformatter_kw : dict-like, optional
            The *x* and *y* axis formatter settings. Passed to
            `~proplot.axistools.Formatter`.
        xrotation, yrotation : None or float, optional
            The rotation for *x* and *y* axis tick labels. Defaults to ``0``
            for normal axes, ``45`` for time axes.

        Note
        ----
        If you plot something with a `numpy`
        `datetime64 <https://docs.scipy.org/doc/numpy/reference/arrays.datetime.html>`__,
        `pandas.Timestamp`, `pandas.DatetimeIndex`, `datetime.date`,
        `datetime.time`, or `datetime.datetime` array as the *x* or *y*-axis
        coordinate, the axis ticks and tick labels will be formatted as dates
        by default.

        See also
        --------
        `BaseAxes.format`, `BaseAxes.smart_update`, `~proplot.axistools.Scale`,
        `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`
        """
        # Background patch basics
        self.patch.set_clip_on(False)
        self.patch.set_zorder(-1)
        self.patch.update(rc.fill({'facecolor': 'axes.facecolor', 'alpha': 'axes.alpha'}))

        # Background hatching (useful where we want to highlight invalid data)
        # TODO: Implement hatching for PolarAxes and map axes?
        hatch = rc['axes.hatch']
        if not self._hatch and hatch: # non-empty
            self._hatch = self.fill_between([0,1], 0, 1, zorder=0, # put in back
                linewidth=0, hatch=hatch, # linewidth affects only patch edge; hatch width controlled by hatch.linewidth
                alpha=rc.get('hatch.alpha'), edgecolor=rc.get('hatch.color'),
                facecolor='none', transform=self.transAxes)
        if self._hatch:
            kw = rc.fill({
                'alpha':     'hatch.alpha',
                'edgecolor': 'hatch.color',
                })
            self._hatch.update(kw)
            self._hatch.set_hatch(hatch)

        # Flexible keyword args, declare defaults
        xmargin       = _default(xmargin, rc['axes.xmargin'])
        ymargin       = _default(ymargin, rc['axes.ymargin'])
        xtickdir      = _default(xtickdir, rc['xtick.direction'])
        ytickdir      = _default(ytickdir, rc['ytick.direction'])
        xtickminor    = _default(xtickminor, rc['xtick.minor.visible'])
        ytickminor    = _default(ytickminor, rc['ytick.minor.visible'])
        xspineloc     = _default(xloc, xspineloc, _rcloc_to_stringloc('x', 'axes.spines'))
        yspineloc     = _default(yloc, yspineloc, _rcloc_to_stringloc('y', 'axes.spines'))
        xformatter    = _default(xticklabels, xformatter) # default is just matplotlib version
        yformatter    = _default(yticklabels, yformatter)
        xlocator      = _default(xticks, xlocator) # default is AutoLocator, no setting
        ylocator      = _default(yticks, ylocator)
        xminorlocator = _default(xminorticks, xminorlocator) # default is AutoMinorLocator, no setting
        yminorlocator = _default(yminorticks, yminorlocator)
        # Grid defaults are more complicated
        # NOTE: rcmod will always change *grid* and *grid.which* at the same time.
        # TODO: Make xgridminor, ygridminor, xtickloc, ytickloc all global rc
        # settings that modify xtick.left, axes.grid.axis, etc., instead of
        # these keyword args. Results in duplicate behavior maybe.
        which = rc['axes.grid.which']
        grid = rc.get('axes.grid', cache=True)
        axis = rc.get('axes.grid.axis', cache=False) # always need this property
        if which is not None or grid is not None: # only if *one* was changed recently!
            # But no matter what we need *both* to figure out proper xgrid, ygrid arguments
            # NOTE: Should we try to make xgridminor, ygridminor part of thing?
            if grid is None:
                grid = rc.get('axes.grid')
            elif which is None:
                which = rc.get('axes.grid.which')
            xgrid      = _default(xgrid, grid and axis in ('x','both') and which in ('major','both'))
            ygrid      = _default(ygrid, grid and axis in ('y','both') and which in ('major','both'))
            xgridminor = _default(xgridminor, grid and axis in ('x','both') and which in ('minor','both'))
            ygridminor = _default(ygridminor, grid and axis in ('y','both') and which in ('minor','both'))

        # Override for weird bug where title doesn't get automatically offset
        # from ticklabels in certain circumstance; check out notebook
        xlabelloc = _default(xlabelloc, xticklabelloc)
        ylabelloc = _default(ylabelloc, yticklabelloc)
        xtickloc = _default(xtickloc, xticklabelloc, _rcloc_to_stringloc('x', 'xtick'))
        ytickloc = _default(ytickloc, yticklabelloc, _rcloc_to_stringloc('y', 'ytick'))
        if xlabelloc=='both':
            xlabelloc = 'bottom'
        if ylabelloc=='both':
            ylabelloc = 'left'
        if xticklabelloc=='both' and xtickloc!='both':
            xticklabelloc = xtickloc
        if yticklabelloc=='both' and ytickloc!='both':
            yticklabelloc = ytickloc
        if xticklabelloc in ('both','top') and (xlabelloc!='top' or not xlabel): # xtickloc *cannot* be 'top', *only* appears for 'both'
            warnings.warn('This keyword combo causes matplotlib bug where title is not offset from tick labels. Try again with xticklabelloc="bottom" or xlabelloc="top". Defaulting to the former.')
            xticklabelloc = xlabelloc = 'bottom'

        # Begin loop
        for (axis, label, color, margin,
            tickloc, spineloc, ticklabelloc, labelloc,
            bounds, grid, gridminor, tickminor, tickminorlocator,
            lim, reverse, scale, locator,
            formatter, tickrange, tickdir, ticklabeldir, rotation,
            scale_kw, label_kw, formatter_kw, locator_kw, minorlocator_kw) in zip(
            (self.xaxis, self.yaxis), (xlabel, ylabel), (xcolor, ycolor), (xmargin, ymargin),
            (xtickloc, ytickloc), (xspineloc, yspineloc), (xticklabelloc, yticklabelloc), (xlabelloc, ylabelloc),
            (xbounds, ybounds), (xgrid, ygrid), (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
            (xlim, ylim), (xreverse, yreverse), (xscale, yscale), (xlocator, ylocator),
            (xformatter, yformatter), (xtickrange, ytickrange), (xtickdir, ytickdir), (xticklabeldir, yticklabeldir), (xrotation, yrotation),
            (xscale_kw, yscale_kw), (xlabel_kw, ylabel_kw), (xformatter_kw, yformatter_kw), (xlocator_kw, ylocator_kw), (xminorlocator_kw, yminorlocator_kw),
            ):
            # Axis label properties
            # Redirect user request to the correct *shared* axes, then
            # to the correct *spanning* axis label.
            name = axis.axis_name
            if label is not None:
                kw = rc.fill({
                    'color':    'axes.edgecolor',
                    'fontname': 'fontname',
                    'fontsize': 'axes.labelsize',
                    'weight':   'axes.labelweight'
                    })
                if axis.get_label_position() == 'top':
                    kw['va'] = 'bottom' # baseline was cramped if no ticklabels present
                if color:
                    kw['color'] = color
                kw.update(label_kw)
                self.figure._axis_label_update(axis, text=label, **kw)

            # Axis scale and limits. These don't have axis-specific setters.
            # If user specified xlocator or ylocator and scale is log, enforce
            # custom formatter; this generally means we want specific tick labels
            # on a log-scale plot, but log formatter overrides this and only shows powers of 10.
            if scale is not None:
                if hasattr(scale, 'name'): # class was passed
                    scale = scale.name
                getattr(self, f'set_{name}scale')(axistools.Scale(scale, **scale_kw))
                if scale in ('log','inverse') and locator is not None and formatter is None:
                    formatter = 'default' # override
            if lim is not None:
                getattr(self, f'set_{name}lim')(lim)
            if reverse:
                # axis.set_inverted(True) # 3.1+, the below is from source code
                lo, hi = axis.get_view_interval()
                axis.set_view_interval(max(lo, hi), min(lo, hi), ignore=True)
            # Detect if time axis
            time = isinstance(axis.converter, mdates.DateConverter) # is this a time axis?

            # Fix spines
            kw = rc.fill({
                'linewidth': 'axes.linewidth',
                'color':     'axes.edgecolor',
                })
            if color is not None:
                kw['color'] = color
            if isinstance(self, mproj.PolarAxes):
                sides = ('inner','polar') if name=='x' else ('start','end')
            else:
                sides = ('bottom','top') if name=='x' else ('left','right')
            spines = [self.spines[s] for s in sides]
            for spine,side in zip(spines,sides):
                # Line properties
                # Override if we're settings spine bounds
                spineloc = getattr(self, name + 'spine_override', spineloc) # optionally override; necessary for twinx/twiny situation
                if bounds is not None and spineloc not in sides:
                    spineloc = sides[0] # by default, should just have spines on edges in this case
                # Eliminate sides
                if spineloc=='neither':
                    spine.set_visible(False)
                elif spineloc=='both':
                    spine.set_visible(True)
                elif spineloc in sides: # make relevant spine visible
                    b = True if side==spineloc else False
                    spine.set_visible(b)
                elif spineloc is not None:
                    # Special spine location
                    # Note special 'spine location' options include 'zero', 'center',
                    # and tuple with (units, location) where units can be axes, data, or outward
                    if side==sides[0]: # move the left/semabottom spine onto the specified location, with set_position
                        spine.set_visible(True)
                        spine.set_position(spineloc)
                    else:
                        spine.set_visible(False)
                # Apply spine bounds
                if bounds is not None and spine.get_visible():
                    spine.set_bounds(*bounds)
                spine.update(kw)
            # Get which spines are visible; needed for setting tick locations
            spines = [side for side,spine in zip(sides,spines) if spine.get_visible()]

            # Rc settings, grid settings, major and minor settings
            # Override is just a "new default", but user can override this
            dict_ = lambda prefix: {
                'grid_color':     prefix + '.color',
                'grid_alpha':     prefix + '.alpha',
                'grid_linewidth': prefix + '.linewidth',
                'grid_linestyle': prefix + '.linestyle',
                }
            override = getattr(self, 'grid_override', None)
            grid = _default(grid, override)
            gridminor = _default(gridminor, override)
            for which,igrid in zip(('major', 'minor'), (grid,gridminor)):
                # Turn grid on
                if igrid is not None:
                    axis.grid(igrid, which=which) # toggle with special global props
                # Grid and tick style
                if which=='major':
                    kw_grid = rc.fill(dict_('grid'))
                else:
                    kw_major = kw_grid
                    kw_grid = rc.fill(dict_('gridminor'))
                    kw_grid.update({key:value for key,value in kw_major.items() if key not in kw_grid})
                # Changed rc settings
                kw = _default(rc[name + 'tick.' + which], {})
                kw.pop('visible', None) # invalid setting
                axis.set_tick_params(which=which, **kw_grid, **kw)

            # Tick and ticklabel properties
            # * Weird issue seems to cause set_tick_params to reset/forget that the grid
            #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
            #   Since gridOn is undocumented feature, don't use it. So calling _format_axes() a second time will remove the lines
            # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
            #   group of left/right or top/bottom
            # * Includes option to draw spines but not draw ticks on that spine, e.g.
            #   on the left/right edges
            # First tick sides
            if not isinstance(self, mproj.PolarAxes):
                kw = {}
                translate = {None: None, 'both': sides, 'neither': (), 'none': ()}
                if bounds is not None and tickloc not in sides:
                    tickloc = sides[0] # override to just one side
                ticklocs = translate.get(tickloc, (tickloc,))
                if ticklocs is not None:
                    kw.update({side: (side in ticklocs) for side in sides})
                kw.update({side: False for side in sides if side not in spines}) # override
                # Tick label sides
                # Will override to make sure only appear where ticks are
                ticklabellocs = translate.get(ticklabelloc, (ticklabelloc,))
                if ticklabellocs is not None:
                    kw.update({f'label{side}': (side in ticklabellocs) for side in sides})
                kw.update({'label' + side: False for side in sides
                    if (side not in spines or (ticklocs is not None and side not in ticklocs))}) # override
                # The axis label side
                if labelloc is None:
                    if ticklocs is not None:
                        options = [side for side in sides if (side in ticklocs and side in spines)]
                        if len(options)==1:
                            labelloc = options[0]
                elif labelloc not in sides:
                    raise ValueError(f'Got labelloc "{labelloc}", valid options are {sides}.')
                # Apply
                axis.set_tick_params(which='both', **kw)
                if labelloc is not None:
                    axis.set_label_position(labelloc)

            # The tick styles
            # First color and size
            kw = rc.fill({
                'labelcolor': 'tick.labelcolor', # new props
                'labelsize':  'tick.labelsize',
                'color':      name + 'tick.color',
                })
            if color:
                kw['color'] = color
                kw['labelcolor'] = color
            # Tick direction and rotation
            if tickdir=='in':
                kw['pad'] = 1 # ticklabels should be much closer
            if ticklabeldir=='in': # put tick labels inside the plot
                tickdir = 'in'
                pad = rc.get(name + 'tick.major.size') + rc.get(name + 'tick.major.pad') + rc.get(name + 'tick.labelsize')
                kw['pad'] = -pad
            if tickdir is not None:
                kw['direction'] = tickdir
            axis.set_tick_params(which='both', **kw)

            # Settings that can't be controlled by set_tick_params
            # Also set rotation here, otherwise get weird alignment
            # See discussion: https://stackoverflow.com/q/11264521/4970632
            kw = rc.fill({'fontname': 'fontname', 'weight':'tick.labelweight'})
            if rotation is not None:
                kw.update({'rotation':rotation})
                if name=='x':
                    kw.update({'ha':'right' if rotation>0 else 'left'})
                setattr(self, f'_{name}rotated', True)
            for t in axis.get_ticklabels():
                t.update(kw)
            # Margins
            if margin is not None:
                self.margins(**{name: margin})

            # Major and minor locator
            # Also automatically detect whether axis is a 'time axis' (i.e.
            # whether user has plotted something with x/y as datetime/date/np.datetime64
            # objects, and matplotlib automatically set the unit converter)
            if locator is not None:
                locator = axistools.Locator(locator, time=time, **locator_kw)
                axis.set_major_locator(locator)
                if isinstance(locator, mticker.IndexLocator):
                    tickminor = False # minor ticks make no sense for 'index' data
            if not tickminor and tickminorlocator is None:
                axis.set_minor_locator(axistools.Locator('null'))
            elif tickminorlocator is not None:
                axis.set_minor_locator(axistools.Locator(tickminorlocator,
                    minor=True, time=time, **minorlocator_kw))

            # Major and minor formatter
            fixedformatfix = False
            if formatter is not None or tickrange is not None:
                # Override
                if tickrange is not None:
                    if formatter not in (None, 'default'):
                        warnings.warn('The tickrange feature requires proplot.ScalarFormatter formatter. Overriding input formatter.')
                    formatter = 'default'
                # Set the formatter
                formatter = axistools.Formatter(formatter, tickrange=tickrange, time=time, **formatter_kw)
                axis.set_major_formatter(formatter)
                if isinstance(formatter, mticker.FixedFormatter): # if locator is MultipleLocator, first tick gets cut off!
                    fixedformatfix = True
            axis.set_minor_formatter(mticker.NullFormatter())

            # Ensure no out-of-bounds ticks! Even set_smart_bounds() fails sometimes.
            # * Using set_bounds also failed, and fancy method overrides did
            #   not work, so instead just turn locators into fixed version
            # * Most locators take no arguments in __call__, and some do not
            #   have tick_values method, so we just call them.
            if fixticks or fixedformatfix or bounds is not None or axis.get_scale()=='cutoff':
                if bounds is None:
                    bounds = getattr(self, f'get_{name}lim')()
                locator = axistools.Locator([x for x in axis.get_major_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_major_locator(locator)
                locator = axistools.Locator([x for x in axis.get_minor_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_minor_locator(locator)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        if (xlim is not None or ylim is not None) and self._inset_parent:
            self.indicate_inset_zoom()
        super().smart_update(**kwargs)

    def dualx(self, offset=0, scale=1, xscale='linear', xlabel=None, **kwargs):
        """As with `~XYAxes.dualy`, but for the *x*-axis. See `~XYAxes.dualy`."""
        parent = self.get_xscale()
        if parent!='linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current "{parent}" scale.')
            self.set_xscale('linear')
        ax = self.twiny()
        if xlabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "xlabel" keyword argument.')
        xscale = axistools.InvertedScaleFactory(xscale)
        ax.format(xscale=xscale, xlabel=xlabel, **kwargs)
        self._dualx_scale = (offset, scale)

    def dualy(self, offset=0, scale=1, yscale='linear', ylabel=None, **kwargs):
        """
        Makes a secondary *y*-axis for denoting equivalent *y*
        coordinates in **alternate units**. Returns nothing.

        Parameters
        ----------
        scale : float, optional
            The constant multiple applied after scaling data with `transform`.
            Defaults to ``1``.
            For example, if your *y*-axis is meters and you
            want kilometers on the other side, use ``scale=1e-3``.
        offset : float, optional
            The constant offset added after multipyling by `scale`.
            Defaults to ``0``.
            For example, if your *y*-axis is Kelvin and you want degrees
            Celsius on the opposite side, use ``offset=-273.15``.
        yscale : str, optional
            The registered scale name used to transform data to the alternate
            units.  Defaults to ``'linear'``.
            For example, if your *y*-axis is wavenumber and you want wavelength on
            the opposite side, use ``yscale='inverse'``. If your-*y* axis
            is height and you want pressure on the opposite side, use
            ``yscale='pressure'`` (and vice versa).
        ylabel : None or str, optional
            The axis label (highly recommended). A warning will be issued if
            this is not supplied.
        **kwargs
            Formats the new axis. Passed to `BaseAxes.format`.

        Note
        ----
        The axis scale `yscale` is used to transform units on the left axis,
        linearly spaced, to units on the right axis. This means the right
        'axis scale' must scale its data with the *inverse* of this transform.
        We make this inverted scale with `~proplot.axistools.InvertedScaleFactory`.
        """
        # Notes:
        # For some reason, when scale is applied, it can change the default
        # formatter. For example, a linear scale will change default formatter
        # to original matplotlib version instead of my custom override. Need
        # to apply it explicitly.
        parent = self.get_yscale()
        if parent!='linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current "{parent}" scale.')
            self.set_yscale('linear')
        ax = self.twinx()
        if ylabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "ylabel" keyword argument.')
        yscale = axistools.InvertedScaleFactory(yscale)
        ax.format(yscale=yscale, ylabel=ylabel, **kwargs)
        self._dualy_scale = (offset, scale)

    def altx(self, *args, **kwargs):
        """Alias (and more intuitive name) for `~XYAxes.twiny`.
        The matplotlib `~matplotlib.axes.Axes.twiny` function
        actually generates two *x*-axes with a shared ("twin") *y*-axis."""
        return self.twiny(*args, **kwargs)

    def alty(self, *args, **kwargs):
        """Alias (and more intuitive name) for `~XYAxes.twinx`.
        The matplotlib `~matplotlib.axes.Axes.twinx` function
        actually generates two *y*-axes with a shared ("twin") *x*-axis."""
        return self.twinx(*args, **kwargs)

    def twinx(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twinx` and intelligently handles
        axis ticks, gridlines, axis tick labels, axis labels, and axis sharing.
        Returns an `XYAxes` instance."""
        # Note: Cannot wrap twinx() because then the axes created will be
        # instantiated from the parent class, which doesn't have format method.
        # Instead, use hidden method _make_twin_axes.
        if self._alty_child:
            raise ValueError('No more than two twin axes!')
        if self._alty_parent:
            raise ValueError('This *is* a twin axes!')
        ax = self._make_twin_axes(sharex=self, projection=self.name)
        # Setup
        self.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_offset_position('right')
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.xaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format called
        self.yspine_override = 'left' # original axis ticks on left
        ax.yspine_override = 'right' # new axis ticks on right
        ax.xspine_override = 'neither'
        ax.grid_override = False
        # Return
        ax._alty_parent = self
        self._alty_child = ax
        return ax

    def twiny(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twiny` and intelligently handles
        axis ticks, gridlines, axis tick labels, axis labels, and axis sharing.
        Returns an `XYAxes` instance."""
        # Note: Cannot wrap twiny() because we want to use our own XYAxes,
        # not the matplotlib Axes. Instead use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py
        if self._altx_child:
            raise ValueError('No more than two twin axes!')
        if self._altx_parent:
            raise ValueError('This *is* a twin axes!')
        ax = self._make_twin_axes(sharey=self, projection=self.name)
        # Setup
        self.xaxis.tick_bottom()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_autoscaley_on(self.get_autoscaley_on())
        ax.yaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format called
        self.xspine_override = 'bottom' # original axis ticks on bottom
        ax.xspine_override = 'top' # new axis ticks on top
        ax.yspine_override = 'neither'
        ax.grid_override = False
        # Return
        ax._altx_parent = self
        self._altx_child = ax
        return ax

    def inset(self, *args, **kwargs):
        """Alias for `~XYAxes.inset_axes`."""
        return self.inset_axes(*args, **kwargs)

    def inset_axes(self, bounds, *, transform=None, zorder=5, zoom=True, zoom_kw={}, **kwargs):
        """Draws an inset `XYAxes` axes. Otherwise, this is a carbon copy
        of the `~matplotlib.axes.Axes.inset_axes` method."""
        # Carbon copy, but use my custom axes
        # Defaults
        if transform is None:
            transform = self.transAxes
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        locator = self._make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = maxes.Axes(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)
        # The following locator lets the axes move if in data coordinates, gets called in ax.apply_aspect()
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        self._inset_children.append(ax)
        ax._inset_parent = self
        # Finally add zoom (NOTE: Requires version >=3.0)
        if zoom:
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def inset_zoom(self, *args, **kwargs):
        """Alias for `~XYAxes.indicate_inset_zoom`."""
        # Just an alias
        return self.indicate_inset_zoom(*args, **kwargs)

    def indicate_inset_zoom(self, alpha=None, linewidth=None, color=None, edgecolor=None, **kwargs):
        """Custom inset zoom indicator that can be *refreshed* as the
        parent axis limits change."""
        # Makes more sense to be defined on the inset axes, since parent
        # could have multiple insets
        parent = self._inset_parent
        alpha = alpha or 1.0
        linewidth = linewidth or rc['axes.linewidth']
        edgecolor = color or edgecolor or rc['axes.edgecolor']
        if not parent:
            raise ValueError(f'{self} is not an inset axes.')
        xlim = self.get_xlim()
        ylim = self.get_ylim()
        rect = [xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0]]
        kwargs.update({'linewidth':linewidth, 'edgecolor':edgecolor, 'alpha':alpha})
        rectpatch, connects = parent.indicate_inset(rect, self, **kwargs)
        # Adopt properties from old one
        if self._zoom:
            rectpatch_old, connects_old = self._zoom
            rectpatch.update_from(rectpatch_old)
            rectpatch_old.set_visible(False)
            for line,line_old in zip(connects,connects_old):
                # Actually want to *preserve* whether line is visible! This
                # is automatically determined!
                visible = line.get_visible()
                line.update_from(line_old)
                line.set_visible(visible)
                line_old.set_visible(False)
        # By default linewidth is only applied to box
        else:
            for line in connects:
                line.set_linewidth(linewidth)
                line.set_color(edgecolor)
                line.set_alpha(alpha)
        self._zoom = (rectpatch, connects)
        return (rectpatch, connects)

    def _make_inset_locator(self, bounds, trans):
        """Helper function, copied from private matplotlib version."""
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
            return bb
        return inset_locator

class EmptyPanel(object):
    """
    Replaces `PanelAxes` when the axes or figure panel does not exist.
    This gives a nicer error message than if we had just ``None``, and
    permits indexing to mimick the behavior of a singleton
    `~proplot.subplots.axes_list`.

    Note
    ----
    `__getattr__` is invoked only when `__getattribute__` fails, i.e.
    when the user requests anything that isn't a builtin method.
    """
    def __bool__(self):
        """Returns False. Provides shorthand way to check whether panel
        attribute is specifically EmptyPanel."""
        return False # it's empty, so this is 'falsey'

    def __len__(self, key):
        """Returns 1. This allows us to treat `EmptyPanel` like an
        `~proplot.subplots.axes_list` of stacked panels."""
        return 1

    def __getitem__(self, key):
        """Returns itself. This allows us to treat `EmptyPanel` like an
        `~proplot.subplots.axes_list` of stacked panels."""
        # See: https://stackoverflow.com/a/26611639/4970632
        if key>0:
            raise IndexError
        return self

    def __getattr__(self, attr, *args):
        """Raises RuntimeError."""
        raise RuntimeError('Panel does not exist.')

class PanelAxes(XYAxes):
    """`~proplot.axes.XYAxes` subclass, adds `~PanelAxes.legend` and
    `~PanelAxes.colorbar` methods that "fill" the entire axes."""
    # Notes:
    # See `this post <https://stackoverflow.com/a/52121237/4970632>`_
    # and `this example <https://stackoverflow.com/q/26236380/4970632>`_.
    name = 'panel'
    """The registered projection name."""
    def __init__(self, *args,
        side=None, share=False, flush=False,
        visible=True, parent=None, **kwargs):
        """
        Parameters
        ----------
        side : {'left', 'right', 'bottom', 'top'}
            The side on which the panel is drawn.
        share : bool, optional
            Whether to share panel *x* and *y* axes with "parent" axes.
            Irrelevant for figure panels, and defaults to ``False``.
        flush : bool, optional
            Whether panel should always be "flush" against its parent
            subplot or not.
        visible : bool, optional
            Whether to make the axes visible or not.
        parent : None or `~matplotlib.axes.Axes`
            The "parent" of the panel. Not relevant for "outer panel"
            axes.
        *args, **kwargs
            Passed to the `XYAxes` initializer.

        See also
        --------
        `~proplot.subplots.subplots`,
        `~proplot.subplots.Figure.add_subplot_and_panels`,
        `XYAxes`, `BaseAxes`
        """
        # Misc props
        # WARNING: Need to set flush property before calling super init!
        self._share = share
        self._side = side
        self._flush = flush # should panel always be flush against subplot?
        self._parent = parent # used when declaring parent
        # Initialize
        super().__init__(*args, **kwargs)
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid panel side "{side}".')
        if not visible:
            self.set_visible(False)

    def legend(self, *args, fill=True, **kwargs):
        """
        Use ``fill=True`` (the default) to "fill the panel" with a legend.
        Use ``fill=False`` to add an inset legend with `matplotlib.axes.Axes.legend`.
        If the former, a centered legend is drawn and the axes spines, ticks,
        etc.  are made invisible.

        See `~matplotlib.axes.Axes.legend` and `legend` for details.
        """
        # Regular old inset legend
        if not fill:
            return super().legend(*args, **kwargs)
        # Hide content
        self._filled = True
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        # Allocate invisible axes for drawing legend; by default try to
        # make handles and stuff flush against the axes edge
        kwdefault = {'borderaxespad': 0}
        if not kwargs.get('frameon', rc.get('legend.frameon')):
            kwdefault['borderpad'] = 0
        kwdefault.update(kwargs)
        kwargs = kwdefault
        # Set location by panel side
        # WARNING: center left and center right also turn off horizontal
        # center alignment, so not an option.
        if 'loc' in kwargs:
            warnings.warn(f'Overriding user input legend property "loc".')
        kwargs['loc'] = {'bottom':'upper center', 'right':'center',
                          'left':'center',  'top':'lower center'}[self._side]
        # For filled axes, call wrapper method directly
        return wrappers.legend_wrapper(self, *args, **kwargs)

    def colorbar(self, *args, fill=True, length=1, **kwargs):
        """
        Use ``fill=True`` (the default) to fill the panel with a colorbar.
        Use ``fill=False`` to add an inset colorbar with `BaseAxes.colorbar`.
        If the former, `length` denotes the fractional extent of the panel
        that is filled.

        See `~matplotlib.figure.Figure.colorbar` and `colorbar`
        for details.
        """
        # Inset 'legend-style' colorbar
        if not fill:
            return super().colorbar(*args, **kwargs)
        # Hide content
        self._filled = True
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        # Draw colorbar with arbitrary length relative to full length of panel
        fig = self.figure
        side = self._side
        subspec = self.get_subplotspec()
        if side=='top': # this is ugly, and hard to implement with title, super title, and stuff
            raise NotImplementedError('Colorbars in upper panels are not allowed.')
        if length!=1:
            if side in ['bottom']:
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=1, ncols=3, wspace=0, #hspace=space,
                        subplot_spec=subspec,
                        width_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subspec = gridspec[1]
            elif side in ['left','right']:
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=3, ncols=1, hspace=0, #wspace=space,
                        subplot_spec=subspec,
                        height_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subspec = gridspec[1]
        # Get properties
        ax = fig.add_subplot(subspec, projection=None)
        if side in ['bottom','top']:
            outside, inside = 'bottom', 'top'
            if side=='top':
                outside, inside = inside, outside
            ticklocation = outside
            orientation  = 'horizontal'
        elif side in ['left','right']:
            outside, inside = 'left', 'right'
            if side=='right':
                outside, inside = inside, outside
            ticklocation = outside
            orientation  = 'vertical'
        # For filled axes, call wrapper method directly
        self._colorbar_child = ax # these are so far unused
        ax._colorbar_parent = self
        kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})
        return wrappers.colorbar_wrapper(ax, *args, **kwargs)

class MapAxes(BaseAxes):
    """Intermediate class, shared by `CartopyAxes` and `BasemapAxes`.
    Disables methods that are inappropriate for map projections and adds
    `MapAxes.smart_update`, so that arguments passed to `~BaseAxes.format` are
    identical for `CartopyAxes` and `BasemapAxes`."""
    def __init__(self, *args, **kwargs): # just to disable docstring inheritence
        """
        See also
        --------
        `~proplot.subplots.subplots`, `BaseAxes`, `CartopyAxes`, `BasemapAxes`
        """
        super().__init__(*args, **kwargs)
        self._is_map = True # needed by wrappers, which can't import this file

    @wrappers._expand_methods_list
    def __getattribute__(self, attr, *args):
        """Disables methods inappropriate for map projections:
        `_map_disabled_methods`."""
        # See: https://stackoverflow.com/a/23126260/4970632
        if attr in wrappers._map_disabled_methods:
            raise NotImplementedError('Invalid plotting function {} for map projection axes.'.format(attr))
        return super().__getattribute__(attr, *args)

    # Note this *actually* just returns some standardized arguments
    # to the CartopyAxes.smart_update and BasemapAxes.smart_update methods; they
    # both jump over this intermediate class and call BaseAxes.smart_update
    def smart_update(self, labels=None, latlabels=None, lonlabels=None,
        latmax=None, lonlim=None, latlim=None, grid=None,
        lonlocator=None, lonlines=None, lonticks=None,
        latlocator=None, latlines=None, latticks=None,
        **kwargs,
        ):
        """
        Format the map projection.

        Parameters
        ----------
        labels : None or bool, optional
            Whether to draw longitude and latitude labels. If ``None``, read
            from `~proplot.rcmod.rc` configuration.
        lonlabels, latlabels
            Whether to label longitudes and latitudes, and on which sides
            of the map. There are four different options:

            1. Boolean ``True``. Indicates left side for latitudes,
               bottom for longitudes.
            2. A string, e.g. ``'lr'`` or ``'bt'``.
            3. A boolean ``(left,right)`` tuple for longitudes,
               ``(bottom,top)`` for latitudes.
            4. A boolean ``(n1,n2,n3,n4)`` tuple as in the
               `~mpl_toolkits.basemap.Basemap.drawmeridians` and
               `~mpl_toolkits.basemap.Basemap.drawparallels` methods.
               The boolean values indicate whether to label gridlines intersecting
               the left, right, top, and bottom sides, respectively.

        latmax : None or float, optional
            Meridian gridlines are cut off poleward of this latitude. If
            ``None``, read from the configuration.
        lonlim, latlim : None or length-2 list of float, optional
            Longitude and latitude limits of projection, applied
            with `~cartopy.mpl.geoaxes.GeoAxes.set_extent`.
        grid : None or bool, optional
            Whether to add gridlines.
        lonlines, latlines, lonticks, latticks
            Aliases for `lonlocator`, `latlocator`.
        lonlocator, latlocator : None or list of float, optional
            List of longitudes and latitudes for drawing gridlines.
        **kwargs
            Passed to `BaseAxes.smart_update`.

        See also
        --------
        `BaseAxes.format`, `BaseAxes.smart_update`, `~proplot.subplots.subplots`, `~proplot.rcmod`
        """
        # Parse alternative keyword args
        # NOTE: If labels keyword args were passed, automatically turn grid on
        grid = _default(grid, rc.get('geogrid'))
        labels = _default(labels, rc.get('geogrid.labels')) or bool(lonlabels or latlabels)
        lonlocator = _default(lonlines, lonticks, lonlocator)
        latlocator = _default(latlines, latticks, latlocator)
        latlocator = _default(latlocator, rc.get('geogrid.latlines')) # gridlines by default
        lonlocator = _default(lonlocator, rc.get('geogrid.lonlines'))

        # Interptet latitude
        latmax = _default(latmax, rc.get('geogrid.latmax'))
        if isinstance(self, CartopyAxes):
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
        else:
            lon_0 = self.m.lonmin
        if lonlocator is not None:
            if not np.iterable(lonlocator):
                lonlocator = utils.arange(lon_0 - 180, lon_0 + 180, lonlocator)
            lonlocator = [*lonlocator]
        if latlocator is not None:
            if not np.iterable(latlocator):
                latlocator = utils.arange(-latmax, latmax, latlocator)
            latlocator = [*latlocator]

        # Length-4 boolean arrays of whether and where to goggle labels
        if lonlabels or latlabels:
            labels = True # toggle them at all?
        ilabels = [lonlabels, latlabels]
        for i,(mode,jlabels) in enumerate(zip(('x', 'y'), tuple(ilabels))):
            if jlabels is False:
                return [0]*4
            if jlabels is None:
                jlabels = 1
                # jlabels = False
                # jlabels = True # will label lons on bottom, lats on left
            if isinstance(jlabels, str):
                string = jlabels
                jlabels = [0]*4
                for idx,char in zip([0,1,2,3],'lrbt'):
                    if char in string:
                        jlabels[idx] = 1
            if isinstance(jlabels, Number): # e.g. *boolean*
                jlabels = np.atleast_1d(jlabels)
            if len(jlabels)==1:
                jlabels = [*jlabels, 0] # default is to label bottom/left
            if len(jlabels)==2:
                if mode=='x':
                    jlabels = [0, 0, *jlabels]
                elif mode=='y':
                    jlabels = [*jlabels, 0, 0]
            elif len(jlabels)!=4:
                raise ValueError(f'Invalid {mode} labels: {jlabels}.')
            ilabels[i] = jlabels
        lonlabels, latlabels = ilabels
        return grid, latmax, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs

class PolarAxes(XYAxes, mproj.PolarAxes):
    """Intermediate class, mixes `~matplotlib.projections.polar.PolarAxes`
    with `MapAxes`."""
    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`, `BaseAxes`, `MapAxes`

        Warning
        -------
        Polar axes have not been tested yet!
        """
        super().__init__(*args, **kwargs)

    def smart_update(self, *args, **kwargs):
        """Calls `BaseAxes.smart_update`."""
        super().smart_update(*args, **kwargs)

    name = 'propolar'
    """The registered projection name."""

# Cartopy takes advantage of documented feature where any class with method
# named _as_mpl_axes can be passed as 'projection' object.
# Feature documented here: https://matplotlib.org/devel/add_new_projection.html
class CartopyAxes(MapAxes, GeoAxes):
    """
    Axes subclass for plotting `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
    projections. Initializes the `cartopy.crs.Projection` instance. Also
    allows for *partial* coverage of azimuthal projections by zooming into
    the full projection, then drawing a circle boundary around some latitude
    away from the center (this is surprisingly difficult to do).
    """
    name = 'cartopy'
    """The registered projection name."""
    # Helper
    _n_bounds = 100 # number of points for drawing circle map boundary
    _proj_np = ('npstere',)
    _proj_sp = ('spstere',)
    _proj_circles = ('laea', 'aeqd', 'stere', 'npstere', 'spstere')
    def __init__(self, *args, map_projection=None, centerlat=90, boundinglat=0, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        centerlat : {90, -90}, optional
            For polar projections, the center latitude of the circle.
        boundinglat : float, optional
            For polar projections, the edge latitude of the circle.
        *args, **kwargs
            Passed to `BaseAxes.__init__`.

        Note
        ----
        The circle stuff for polar projection was developed from `this example
        <https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`_.

        See also
        --------
        `~proplot.proj`, `BaseAxes`, `MapAxes`, `~proplot.subplots.subplots`
        """
        # Dependencies
        import cartopy.crs as ccrs # verify package is available

        # Do the GeoAxes initialization steps manually (there are very few)
        # If _hold is set False or None, cartopy will call cla() on axes before
        # plotting stuff, which will wipe out row and column labels even though
        # they appear to stick around; maybe the artist persists but is no longer
        # associated with the axes. Does not matter whether attribute is hidden.
        # self._hold = None
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError('You must initialize CartopyAxes with map_projection=<cartopy.crs.Projection>.')
        self._cartopy_gl = None # gridliner
        self.projection = map_projection # attribute used extensively by GeoAxes methods, and by builtin one

        # Call BaseAxes
        super().__init__(*args, map_projection=map_projection, **kwargs)
        # Apply circle boundary
        name = map_projection.proj4_params['proj']
        if name not in self._proj_circles:
            self.set_global() # see: https://stackoverflow.com/a/48956844/4970632
        else:
            if name in self._proj_np:
                centerlat = 90
            elif name in self._proj_sp:
                centerlat = -90
            eps = 1e-3 # had bug with full -180, 180 range when lon_0 was not 0
            center = self.projection.proj4_params['lon_0']
            self.set_extent([center - 180 + eps, center + 180 - eps, boundinglat, centerlat], ccrs.PlateCarree()) # use platecarree transform
            self.set_boundary(projs.Circle(self._n_bounds), transform=self.transAxes)

    def __getattribute__(self, attr, *args):
        """Wraps methods when they are requested by the user. See
        `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`,
        `~proplot.wrappers.cartopy_transform`, `~proplot.wrappers.check_centers`,
        `~proplot.wrappers.check_edges`, and `~proplot.wrappers.cartopy_gridfix`
        for more info."""
        obj = super().__getattribute__(attr, *args)
        # Color wrappers
        if attr in wrappers._cmap_methods:
            obj = wrappers._cmap_wrapper(self, obj)
        elif attr in wrappers._cycle_methods:
            obj = wrappers._cycle_wrapper(self, obj)
        # Plotting wrappers
        if attr=='plot':
            obj = wrappers._plot_wrapper(self, obj)
        elif attr=='scatter':
            obj = wrappers._scatter_wrapper(self, obj)
        elif attr in wrappers._edges_methods or attr in wrappers._centers_methods:
            obj = wrappers._cartopy_gridfix(self, obj)
            if attr in wrappers._edges_methods:
                obj = wrappers._check_edges(self, obj)
            else:
                obj = wrappers._check_centers(self, obj)
        # Standardized input
        if attr in wrappers._transform_methods:
            obj = wrappers._cartopy_transform(self, obj)
        elif attr in wrappers._crs_methods:
            obj = wrappers._cartopy_crs(self, obj)
        return obj

    def smart_update(self, grid=None, **kwargs):
        # Documentation inherited from MapAxes
        # Dependencies
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs

        # Parse flexible input
        grid, _, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs = \
                super().smart_update(**kwargs)

        # Projection extent
        # NOTE: They may add this as part of set_xlim and set_ylim in the
        # near future; see: https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638
        # WARNING: The set_extent method tries to set a *rectangle* between
        # the *4* (x,y) coordinate pairs (each corner), so something like
        # (-180,180,-90,90) will result in *line*, causing error!
        if lonlim is not None or latlim is not None:
            lonlim = lonlim or [None, None]
            latlim = latlim or [None, None]
            lonlim, latlim = [*lonlim], [*latlim]
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
            if lonlim[0] is None:
                lonlim[0] = lon_0 - 180
            if lonlim[1] is None:
                lonlim[1] = lon_0 + 180
            eps = 1e-3 # had bug with full -180, 180 range when lon_0 was not 0
            lonlim[0] += eps
            if latlim[0] is None:
                latlim[0] = -90
            if latlim[1] is None:
                latlim[1] = 90
            self.set_extent([*lonlim, *latlim], ccrs.PlateCarree())

        # Draw gridlines
        # Make old one invisible
        if self._cartopy_gl is not None:
            gl = self._cartopy_gl
            gl.xlines         = False
            gl.xlabels_top    = False
            gl.xlabels_bottom = False
            gl.ylines         = False
            gl.ylabels_left   = False
            gl.ylabels_right  = False
        # Make new gridliner
        # WARNING: For some reason very weird side effects happen if you try
        # to call gridlines twice on same axes. So impossible to do 'major'
        # and 'minor' gridlines.
        if grid:
            kw = rc.fill({
                'alpha':     'geogrid.alpha',
                'color':     'geogrid.color',
                'linewidth': 'geogrid.linewidth',
                'linestyle': 'geogrid.linestyle',
                }, cache=False)
            labels = labels and isinstance(self.projection, (ccrs.Mercator, ccrs.PlateCarree))
            gl = self.gridlines(draw_labels=labels, zorder=100, **kw)
            self._cartopy_gl = gl

            # Grid locations
            eps = 1e-3
            if latlocator[0]==-90:
                latlocator[0] += eps
            if latlocator[-1]==90:
                latlocator[-1] -= eps
            gl.ylocator = mticker.FixedLocator(latlocator)
            gl.xlocator = mticker.FixedLocator(lonlocator)

            # Grid labels
            if labels:
                from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
                self._gridliner_on = True
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabels_bottom, gl.xlabels_top = lonlabels[2:]
                gl.ylabels_left, gl.ylabels_right = latlabels[:2]
            else:
                self._gridliner_on = False

        # Geographic features
        # WARNING: Seems cartopy features can't be updated!
        # See: https://scitools.org.uk/cartopy/docs/v0.14/_modules/cartopy/feature.html#Feature
        # Change the _kwargs property also does *nothing*
        # WARNING: Changing linewidth is evidently impossible with cfeature. Bug?
        # See: https://stackoverflow.com/questions/43671240/changing-line-width-of-cartopy-borders
        # NOTE: The natural_earth_shp method is deprecated, use add_feature instead.
        # See: https://cartopy-pelson.readthedocs.io/en/readthedocs/whats_new.html
        # NOTE: The e.g. cfeature.COASTLINE features are just for convenience,
        # hi res versions. Use cfeature.COASTLINE.name to see how it can be looked
        # up with NaturalEarthFeature.
        reso = rc.reso
        if reso not in ('lo','med','hi'):
            raise ValueError(f'Invalid resolution {reso}.')
        reso = {
            'lo':  '110m',
            'med': '50m',
            'hi':  '10m',
            }.get(reso)
        features = {
            'land':         ('physical', 'land'),
            'ocean':        ('physical', 'ocean'),
            'lakes':        ('physical', 'lakes'),
            'coast':        ('physical', 'coastline'),
            'rivers':       ('physical', 'rivers_lake_centerlines'),
            'borders':      ('cultural', 'admin_0_boundary_lines_land'),
            'innerborders': ('cultural', 'admin_1_states_provinces_lakes'),
            }
        for name, args in features.items():
            # Get feature
            if not rc.get(name): # toggled
                continue
            if getattr(self, f'_{name}', None): # already drawn
                continue
            feat = cfeature.NaturalEarthFeature(*args, reso)
            # Customize
            # For 'lines', need to specify edgecolor and facecolor individually
            # See: https://github.com/SciTools/cartopy/issues/803
            kw = {**rc.get(name, nodict=False)}
            if name in ('coast', 'rivers', 'borders', 'innerborders'):
                kw['edgecolor'] = kw.pop('color')
                kw['facecolor'] = 'none'
            else:
                kw['linewidth'] = 0
            self.add_feature(feat, **kw)
            setattr(self, f'_{name}', feat)

        # Update patch
        kw = rc.fill({
            'facecolor': 'axes.facecolor'
            }, cache=False)
        self.background_patch.update(kw)
        kw = rc.fill({
            'edgecolor': 'axes.edgecolor',
            'linewidth': 'axes.linewidth'
            }, cache=False)
        self.outline_patch.update(kw)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        BaseAxes.smart_update(self, **kwargs)

def _ls_translate(obj, style):
    """Make basemap gridlines look like cartopy lines using the `dashes`
    tuple."""
    # Notes:
    # * For some reason basemap gridlines look different from cartopy ones.
    #   Have absolutely **no idea** why. Cartopy seems to do something weird because
    #   there is no ``_dashSeq`` attribute on the lines, and the line styles
    #   are always "officially" ``'-'``, even though they are dashed.
    # * Dots actually look better for cartopy so we mimick them. After testing,
    #   below works *really* well... maybe it is hardcoded in cartopy somewhere?
    # See `this reference <https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html>`_.
    dashes = [None, None]
    if style==':':
        dashes[0] = 0.05
        dashes[1] = 2.5
    elif style=='--':
        dashes[0] = 2.5
        dashes[1] = 2.5
    elif style!='-':
        raise ValueError(f'Invalid line style {style}.')
    return dashes

class BasemapAxes(MapAxes):
    """
    Axes subclass for plotting `basemap <https://matplotlib.org/basemap/>`_
    projections. The `~mpl_toolkits.basemap.Basemap` projection instance is added as
    the `m` attribute, but this is all abstracted away -- you can use
    `~matplotlib.axes.Axes` methods like `~matplotlib.axes.Axes.plot` and
    `~matplotlib.axes.Axes.contour` with
    your raw longitude-latitude data.
    """
    name = 'basemap'
    """The registered projection name."""
    # Note non-rectangular projections; for rectnagular ones, axes spines are
    # used as boundaries, but for these, have different boundary.
    _proj_non_rectangular = (
            # Always non-rectangular
            'ortho', 'geos', 'nsper',
            'moll', 'hammer', 'robin',
            'eck4', 'kav7', 'mbtfpq', # last one is McBryde-Thomas flat polar quartic
            'sinu', 'vandg', # last one is van der Grinten
            # Only non-rectangular if we pass 'round' kwarg
            # This is done by default, currently no way to change it
            'npstere', 'spstere', 'nplaea',
            'splaea', 'npaeqd', 'spaeqd',
            )
    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        **kwargs
            Passed to `BaseAxes.__init__`.

        See also
        --------
        `~proplot.proj`, `BaseAxes`, `MapAxes`, `~proplot.subplots.subplots`
        """
        # Some notes
        # * Must set boundary before-hand, otherwise the set_axes_limits method called
        #   by mcontourf/mpcolormesh/etc draws two mapboundary Patch objects called "limb1" and
        #   "limb2" automatically: one for fill and the other for the edges
        # * Then, since the patch object in _mapboundarydrawn is only the fill-version, calling
        #   drawmapboundary again will replace only *that one*, but the original visible edges
        #   are still drawn -- so e.g. you can't change the color
        # * If you instead call drawmapboundary right away, _mapboundarydrawn will contain
        #   both the edges and the fill; so calling it again will replace *both*
        import mpl_toolkits.basemap as mbasemap # verify package is available
        if not isinstance(map_projection, mbasemap.Basemap):
            raise ValueError('You must initialize BasemapAxes with map_projection=(basemap.Basemap instance).')
        self.m = map_projection
        self.boundary = None
        self._mapboundarydrawn = None
        self._recurred = False # use this so we can override plotting methods
        self._parallels = None
        self._meridians = None

        # Initialize
        super().__init__(*args, **kwargs)

    # Basemap overrides
    # WARNING: Never ever try to just make blanket methods on the Basemap
    # instance accessible from axes instance! Can of worms and had bunch of
    # weird errors! Just pick the ones you think user will want to use.
    # NOTE: Need to again use cmap and cycle wrapper here, because
    # Basemap internally bypasses my BaseAxes superclass.
    def __getattribute__(self, attr, *args):
        """
        Wraps methods when they are requested by the user. See
        `~proplot.wrappers.basemap_latlon`, `~proplot.wrappers.check_centers`,
        `~proplot.wrappers.check_edges`, and `~proplot.wrappers.basemap_gridfix` for more info.

        Wraps all plotting methods with ``_wrapper_m_call`` and
        ``_wrapper_m_norecurse``, which (former) calls methods on the
        `~mpl_toolkits.basemap.Basemap` instance and (latter) prevents
        recursion issues arising from internal calls to axes methods by
        the `~mpl_toolkits.basemap.Basemap` instance.
        """
        obj = super().__getattribute__(attr, *args)
        if attr in wrappers._latlon_methods or attr in wrappers._edges_methods \
                or attr in wrappers._centers_methods:
            # Call Basemap method at bottom level
            obj = wrappers._wrapper_m_call(self, obj)
            # Color wrappers
            if attr in wrappers._cmap_methods:
                obj = wrappers._cmap_wrapper(self, obj)
            elif attr in wrappers._cycle_methods:
                obj = wrappers._cycle_wrapper(self, obj)
            # Plotting wrappers
            if attr=='plot':
                obj = wrappers._plot_wrapper(self, obj)
            elif attr=='scatter':
                obj = wrappers._scatter_wrapper(self, obj)
            elif attr in wrappers._edges_methods or attr in wrappers._centers_methods:
                obj = wrappers._basemap_gridfix(self, obj)
                if attr in wrappers._edges_methods:
                    obj = wrappers._check_edges(self, obj)
                else:
                    obj = wrappers._check_centers(self, obj)
            # Standardized input wrappers
            if attr in wrappers._latlon_methods:
                obj = wrappers._basemap_latlon(self, obj)
            # Recursion fix at top level
            obj = wrappers._wrapper_m_norecurse(self, obj)
        return obj

    def smart_update(self, grid=None, **kwargs):
        # Documentation inherited from MapAxes
        # Parse flexible input
        grid, latmax, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs = \
                super().smart_update(**kwargs)
        if lonlim is not None or latlim is not None:
            warnings.warn('You cannot "zoom into" a basemap projection after creating it. Pass a proj_kw dictionary in your call to subplots, with any of the following basemap keywords: llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, llcrnrx, llcrnry, urcrnrx, urcrnry, width, or height.')

        # Map boundary
        # * First have to *manually replace* the old boundary by just deleting
        #   the original one
        # * If boundary is drawn successfully should be able to call
        #   self.m._mapboundarydrawn.set_visible(False) and edges/fill color disappear
        # * For now will enforce that map plots *always* have background whereas
        #   axes plots can have transparent background
        kw_face = rc.fill({
            'facecolor': 'axes.facecolor'
            }, cache=False)
        kw_edge = rc.fill({
            'linewidth': 'axes.linewidth',
            'edgecolor': 'axes.edgecolor'
            }, cache=False)
        self.axesPatch = self.patch # bugfix or something
        if self.m.projection in self._proj_non_rectangular:
            self.patch.set_alpha(0) # make patch invisible
            if not self.m._mapboundarydrawn:
                p = self.m.drawmapboundary(ax=self) # set fill_color to 'none' to make transparent
            else:
                p = self.m._mapboundarydrawn
            p.update({**kw_face, **kw_edge})
            p.set_rasterized(False) # not sure about this; might be rasterized
            p.set_clip_on(False)    # so edges of *line* denoting boundary aren't cut off
            self.boundary = p       # not sure why this one
        else:
            self.patch.update({**kw_face, 'edgecolor':'none'})
            for spine in self.spines.values():
                spine.update(kw_edge)

        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        if grid:
            lkw = rc.fill({
                'alpha':     'geogrid.alpha',
                'color':     'geogrid.color',
                'linewidth': 'geogrid.linewidth',
                'linestyle': 'geogrid.linestyle',
                }, cache=False)
            tkw = rc.fill({
                'color':    'geogrid.color',
                'fontsize': 'geogrid.labelsize',
                }, cache=False)
            # Remove old ones
            if self._parallels:
                for pi in self._parallels.values():
                    for obj in [i for j in pi for i in j]: # magic
                        obj.set_visible(False)
            # Change from left/right/bottom/top to left/right/top/bottom
            if labels:
                latlabels[2:] = latlabels[2:][::-1]
            else:
                latlabels = 4*[0]
            p = self.m.drawparallels(latlocator, latmax=latmax, labels=latlabels, ax=self)
            for pi in p.values(): # returns dict, where each one is tuple
                # Tried passing clip_on to the below, but it does nothing; must set
                # for lines created after the fact
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
                        obj.set_dashes(_ls_translate(obj, lkw['linestyle']))
            self._parallels = p

            # Longitudes
            # Remove old ones
            if self._meridians:
                for pi in self._meridians.values():
                    for obj in [i for j in pi for i in j]: # magic
                        obj.set_visible(False)
            # Draw new ones
            if labels:
                lonlabels[2:] = lonlabels[2:][::-1]
            else:
                lonlabels = 4*[0]
            p = self.m.drawmeridians(lonlocator, latmax=latmax, labels=lonlabels, ax=self)
            for pi in p.values():
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
                        obj.set_dashes(_ls_translate(obj, lkw['linestyle']))
            self._meridians = p

        # Geography
        # NOTE: Also notable are drawcounties, blumarble, drawlsmask,
        # shadedrelief, and etopo methods.
        features = {
            'land':      'fillcontinents',
            'coast':     'drawcoastlines',
            'rivers':    'drawrivers',
            'borders':   'drawcountries',
            'innerborders': 'drawstates',
            }
        for name, method in features.items():
            if not rc.get(name): # toggled
                continue
            if getattr(self, f'_{name}', None): # already drawn
                continue
            feat = getattr(self.m, method)(ax=self, **rc.get(name, nodict=False))
            setattr(self, f'_{name}', feat)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        BaseAxes.smart_update(self, **kwargs)

# Register the projections
mproj.register_projection(BaseAxes)
mproj.register_projection(XYAxes)
mproj.register_projection(PanelAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)

