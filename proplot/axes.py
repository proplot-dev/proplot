#!/usr/bin/env python3
"""
This page documents the axes subclasses returned by
`~proplot.subplots.subplots` and their various method wrappers. You should
start with the documentation on the following methods.

* `Axes.format`
* `Axes.context`
* `CartesianAxes.format`
* `ProjectionAxes.format`

`Axes.format` and `Axes.context` are both called by
`CartesianAxes.format` and `ProjectionAxes.format`. ``format`` is your
**one-stop-shop for changing axes settings** like *x* and *y* axis limits,
axis labels, tick locations, tick labels grid lines, axis scales, titles,
a-b-c labelling, adding geographic features, and much more.

.. raw:: html

   <h1>Developer notes</h1>

Axes method wrappers are documented in the "functions" table. The wrappers are
dynamically applied within the `~proplot.axes.Axes.__getattribute__` methods
on `~proplot.axes.Axes` and its subclasses. But why doesn't ProPlot just
use decorators? *Brevity*. For example, `~proplot.wrappers.cmap_wrapper`
overrides *a dozen* different methods. This lets ProPlot override these
methods in *one* line, instead of 50 lines. To see which methods are overriden,
the user can simply check the documentation.

It should be noted that dynamically wrapping every time the user requests a
method is slower than "decoration", which just wraps the method when the class
is declared. But this was not found to significantly affect performance. And
anyway, `Premature Optimization is the Root of All Evil
<http://wiki.c2.com/?PrematureOptimization>`__.
"""
# WARNING: Wanted to bulk wrap methods using __new__ on a *metaclass*, since
# this would just wrap method *once* and not every single time user accesses
# object. More elegant, but __new__ *does not receive inherited methods* (that
# comes later down the line), so we can't wrap them. Anyway overriding
# __getattribute__ is fine, and premature optimiztaion is root of all evil!
import re
import numpy as np
import warnings
from numbers import Number, Integral
import matplotlib.projections as mproj
import matplotlib.axes as maxes
import matplotlib.dates as mdates
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections
from .rctools import rc, RC_NAMES_NODOTS
from . import utils, projs, axistools, wrappers
from .utils import _notNone, units
__all__ = [
    'Axes',
    'BasemapAxes',
    'CartesianAxes',
    'CartopyAxes',
    'PolarAxes', 'ProjectionAxes',
    ]

# Translator for inset colorbars and legends
SIDE_TRANSLATE = {
    'l':'left',
    'r':'right',
    'b':'bottom',
    't':'top',
    }
LOC_TRANSLATE = {
    None:None,
    'l':'left',
    'r':'right',
    'b':'bottom',
    't':'top',
    'c':'center',
    'i':'best',
    'inset':'best',
    0:'best',
    1:'upper right',
    2:'upper left',
    3:'lower left',
    4:'lower right',
    5:'center left',
    6:'center right',
    7:'lower center',
    8:'upper center',
    9:'center',
    'ur':'upper right',
    'ul':'upper left',
    'll':'lower left',
    'lr':'lower right',
    'cr':'center right',
    'cl':'center left',
    'uc':'upper center',
    'lc':'lower center',
    }

# Helper function
ABC_STRING = 'abcdefghijklmnopqrstuvwxyz'
def _abc(i):
    """Function for a-b-c labeling, returns a...z...aa...zz...aaa...zzz."""
    if i < 26:
        return ABC_STRING[i]
    else:
        return _abc(i - 26) + ABC_STRING[i % 26] # sexy sexy recursion

# Import mapping toolbox
try:
    from cartopy.mpl.geoaxes import GeoAxes
except ModuleNotFoundError:
    GeoAxes = object

#-----------------------------------------------------------------------------#
# Generalized custom axes class
#-----------------------------------------------------------------------------#
class Axes(maxes.Axes):
    """Lowest-level axes subclass. Handles titles and axis
    sharing. Adds several new methods and overrides existing ones."""
    def __init__(self, *args, number=None,
        sharex=None, sharey=None, sharex_level=0, sharey_level=0,
        spanx=None, spany=None, alignx=None, aligny=None,
        **kwargs):
        """
        Parameters
        ----------
        number : int
            The subplot number, used for a-b-c labelling (see
            `~Axes.format`).
        sharex, sharey : `Axes`, optional
            Axes to use for *x* and *y* axis sharing.
        sharex_level, sharey_level : {3, 2, 1, 0}, optional
            The "axis sharing level" for the *x* axis, *y* axis, or both
            axes. See `~proplot.subplots.subplots` for details.
        spanx, spany : bool, optional
            Boolean toggle for whether spanning labels are enabled for the
            *x* and *y* axes. See `~proplot.subplots.subplots` for details.
        alignx, aligny : bool, optional
            Boolean toggle for whether aligned axis labels are enabled for the
            *x* and *y* axes. See `~proplot.subplots.subplots` for details.

        See also
        --------
        `~proplot.subplots.subplots`, `CartesianAxes`, `ProjectionAxes`
        """
        # Call parent
        super().__init__(*args, **kwargs)
        # Properties
        self._number = number # for abc numbering
        self._abc_loc = None
        self._abc_text = None
        self._titles_dict = {} # dictionary of title text objects and their locations
        self._title_loc = None # location of main title
        self._title_pad = rc.get('axes.titlepad') # so we can copy to top panel
        self._title_above_panel = True # TODO: add rc prop?
        # Children and related properties
        self._bpanels = []
        self._tpanels = []
        self._lpanels = []
        self._rpanels = []
        self._tight_bbox = None # bounding boxes are saved
        self._panel_side = None
        self._panel_parent = None
        self._panel_filled = False # True when panels "filled" with colorbar/legend
        self._inset_parent = None
        self._inset_zoom = False
        self._inset_zoom_data = None
        self._alty_child = None
        self._altx_child = None
        self._alty_parent = None
        self._altx_parent = None
        self._auto_colorbar = {} # stores handles and kwargs for auto colorbar
        self._auto_legend = {}
        # Axis sharing, new text attributes, custom formatting
        self._spanx = spanx # boolean toggles, whether we want to span axes labels
        self._spany = spany
        self._alignx = alignx
        self._aligny = aligny
        self._sharex_level = sharex_level
        self._sharey_level = sharey_level
        self._sharex_setup(sharex, sharex_level)
        self._sharey_setup(sharey, sharey_level)
        coltransform = mtransforms.blended_transform_factory(self.transAxes, self.figure.transFigure)
        rowtransform = mtransforms.blended_transform_factory(self.figure.transFigure, self.transAxes)
        self._llabel   = self.text(0.05, 0.5, '', va='center', ha='right', transform=rowtransform)
        self._rlabel  = self.text(0.95, 0.5, '', va='center', ha='left', transform=rowtransform)
        self._blabel = self.text(0.5, 0.05, '', va='top', ha='center', transform=coltransform)
        self._tlabel    = self.text(0.5, 0.95, '', va='bottom', ha='center', transform=coltransform) # reasonable starting point
        self.format(mode=1) # mode == 1 applies the rcExtraParams

    @wrappers._expand_methods_list
    def __getattribute__(self, attr, *args):
        """Disables the redundant methods `DISABLED_METHODS` with useful
        error messages, and applies the `~proplot.wrappers.text_wrapper`
        wrapper."""
        obj = object.__getattribute__(self, attr, *args)
        for message,attrs in wrappers.DISABLED_METHODS.items():
            if attr in attrs:
                raise AttributeError(message.format(attr))
        if attr == 'text':
            obj = wrappers._text_wrapper(self, obj)
        return obj

    def _draw_auto_legends_colorbars(self):
        """Generate automatic legends and colorbars. Wrapper funcs
        let user add handles to location lists with successive calls to
        make successive calls to plotting commands."""
        for loc,(handles,kwargs) in self._auto_colorbar.items():
            self.colorbar(handles, **kwargs)
        for loc,(handles,kwargs) in self._auto_legend.items():
            self.legend(handles, **kwargs)
        self._auto_legend = {}
        self._auto_colorbar = {}

    def _get_side_axes(self, side):
        """Returns groups of axes in row or column or the single group in the
        same row or column as this axes."""
        s = side[0]
        if s not in 'lrbt':
            raise ValueError(f'Invalid side {side!r}.')
        x = ('x' if s in 'lr' else 'y')
        idx = (0 if s in 'lt' else 1) # which side of range to test
        coord = self._range_gridspec(x)[idx] # side for a particular axes
        return [ax for ax in self.figure._axes_main
                if ax._range_gridspec(x)[idx] == coord]

    def _get_title_props(self, abc=False, loc=None):
        """Returns standardized location name, position keyword arguments, and
        setting keyword arguments for the relevant title or a-b-c label at location
        `loc`."""
        # Props
        # NOTE: Sometimes we load all properties from rc object, sometimes
        # just changed ones. This is important if e.g. user calls in two
        # lines ax.format(titleweight='bold') then ax.format(title='text'),
        # don't want to override custom setting with rc default setting.
        props = lambda cache: rc.fill({
            'fontsize':   f'{prefix}.fontsize',
            'weight':     f'{prefix}.weight',
            'color':      f'{prefix}.color',
            'border':     f'{prefix}.border',
            'linewidth':  f'{prefix}.linewidth',
            'fontfamily': 'font.family',
            }, cache=cache)

        # Location string and position coordinates
        cache = True
        prefix = 'abc' if abc else 'title'
        loc = _notNone(loc, rc[f'{prefix}.loc'])
        iloc = getattr(self, '_' + ('abc' if abc else 'title') + '_loc') # old
        if loc is None:
            loc = iloc
        elif iloc is not None and loc != iloc:
            cache = False

        # Above axes
        loc = LOC_TRANSLATE.get(loc, loc)
        if loc in ('top','bottom'):
            raise ValueError(f'Invalid title location {loc!r}.')
        elif loc in ('left','right','center'):
            kw = props(cache)
            kw.pop('border', None) # no border for titles outside axes
            kw.pop('linewidth', None)
            if loc == 'center':
                obj = self.title
            else:
                obj = getattr(self, '_' + loc + '_title')
        # Inside axes
        elif loc in self._titles_dict:
            kw = props(cache)
            obj = self._titles_dict[loc]
        else:
            kw = props(False)
            width, height = self.get_size_inches()
            if loc in ('upper center','lower center'):
                x, ha = 0.5, 'center'
            elif loc in ('upper left','lower left'):
                xpad = rc.get('axes.titlepad')/(72*width)
                x, ha = 1.5*xpad, 'left'
            elif loc in ('upper right','lower right'):
                xpad = rc.get('axes.titlepad')/(72*width)
                x, ha = 1 - 1.5*xpad, 'right'
            else:
                raise ValueError(f'Invalid title or abc "loc" {loc}.')
            if loc in ('upper left','upper right','upper center'):
                ypad = rc.get('axes.titlepad')/(72*height)
                y, va = 1 - 1.5*ypad, 'top'
            elif loc in ('lower left','lower right','lower center'):
                ypad = rc.get('axes.titlepad')/(72*height)
                y, va = 1.5*ypad, 'bottom'
            obj = self.text(x, y, '', ha=ha, va=va, transform=self.transAxes)
            obj.set_transform(self.transAxes)
        return loc, obj, kw

    @staticmethod
    def _loc_translate(loc, **kwargs):
        """Translates location string `loc` into a standardized form."""
        if loc is True:
            loc = None
        elif isinstance(loc, (str, Integral)):
            loc = LOC_TRANSLATE.get(loc, loc) # may still be invalid
        return loc

    def _make_inset_locator(self, bounds, trans):
        """Helper function, copied from private matplotlib version."""
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
            return bb
        return inset_locator

    def _range_gridspec(self, x):
        """Gets the column or row range for the axes. Use `topmost` to get
        properties for the main gridspec grid."""
        subplotspec = self.get_subplotspec()
        if x == 'x':
            _, _, _, _, col1, col2 = subplotspec.get_active_rows_columns()
            return col1, col2
        else:
            _, _, row1, row2, _, _ = subplotspec.get_active_rows_columns()
            return row1, row2

    def _range_tightbbox(self, x):
        """Gets span of tight bounding box, including twin axes and panels
        which are not considered real children and so aren't ordinarily included in
        the tight bounding box calc. `~proplot.axes.Axes.get_tightbbox` caches
        tight bounding boxes when `~Figure.get_tightbbox` is called."""
        # TODO: Better resting for axes visibility
        if x == 'x':
            return self._tight_bbox.xmin, self._tight_bbox.xmax
        else:
            return self._tight_bbox.ymin, self._tight_bbox.ymax

    def _reassign_suplabel(self, side):
        """Re-assigns the column and row labels to panel axes, if they exist.
        This is called by `~proplot.subplots.Figure._align_suplabel`."""
        # Place column and row labels on panels instead of axes -- works when
        # this is called on the main axes *or* on the relevant panel itself
        # TODO: Mixed figure panels with super labels? How does that work?
        s = side[0]
        side = SIDE_TRANSLATE[s]
        if s == self._panel_side:
            ax = self._panel_parent
        else:
            ax = self
        paxs = getattr(ax, '_' + s + 'panels')
        if not paxs:
            return ax
        idx = (0 if s in 'lt' else -1)
        pax = paxs[idx]
        kw = {}
        obj = getattr(ax, '_' + s + 'label')
        for key in ('color', 'fontproperties'): # TODO: add to this?
            kw[key] = getattr(obj, 'get_' + key)()
        pobj = getattr(pax, '_' + s + 'label')
        pobj.update(kw)
        text = obj.get_text()
        if text:
            obj.set_text('')
            pobj.set_text(text)
        return pax

    def _reassign_title(self):
        """Re-assigns title to the first upper panel if present. We cannot
        simply add upper panel as child axes, because then title will be offset
        but still belong to main axes, which messes up tight bounding box."""
        # Reassign title from main axes to top panel -- works when this is
        # called on the main axes *or* on the top panel itself. This is
        # critical for bounding box calcs; not always clear whether draw() and
        # get_tightbbox() are called on the main axes or panel first
        if self._panel_side == 'top':
            ax, taxs = self._panel_parent, [self]
        else:
            ax, taxs = self, self._tpanels
        if not taxs or not ax._title_above_panel:
            tax = ax
        else:
            tax = taxs[0]
            tax._title_pad = ax._title_pad
            for loc,obj in ax._titles_dict.items():
                if not obj.get_text() or loc not in ('left','center','right'):
                    continue
                kw = {}
                loc, tobj, _ = tax._get_title_props(loc=loc)
                for key in ('text', 'color', 'fontproperties'): # TODO: add to this?
                    kw[key] = getattr(obj, 'get_' + key)()
                tobj.update(kw)
                tax._titles_dict[loc] = tobj
                obj.set_text('')

        # Push title above tick marks -- this is known matplotlib problem,
        # but especially annoying with top panels!
        # TODO: Make sure this is robust. Seems 'default' is returned usually
        # when tick label sides is actually *both*. Also makes sure axis is
        # visible; if not, this is a filled colorbar or legend, no padding needed
        pad = 0
        pos = tax.xaxis.get_ticks_position()
        labs = tax.xaxis.get_ticklabels()
        if pos == 'default' or (pos == 'top' and not len(labs)) or (
            pos == 'unknown' and tax._panel_side == 'top'
            and not len(labs) and tax.xaxis.get_visible()):
            pad = tax.xaxis.get_tick_padding()
        tax._set_title_offset_trans(self._title_pad + pad)

    def _share_short_axis(self, share, side, level):
        """When sharing main subplots, shares the short axes of their side
        panels."""
        # TODO: Re-calculate share settings at draw time!
        if self._panel_side: # not None
            return
        s = side[0]
        if s not in 'lrbt':
            raise ValueError(f'Invalid side {side!r}.')
        paxs1 = getattr(self, '_' + s + 'panels') # calling this means, share properties on this axes with input 'share' axes
        paxs2 = getattr(share, '_' + s + 'panels')
        if (not all(not pax._panel_filled for pax in paxs1)
            or not all(not pax._panel_filled for pax in paxs2)):
            return
        if len(paxs1) != len(paxs2):
            raise AttributeError('Sync error. Different number of stacked panels along axes on like column/row of figure.')
        axis = 'x' if s in 'lr' else 'y'
        for pax1,pax2 in zip(paxs1,paxs2):
            getattr(pax1, '_share' + axis + '_setup')(pax2, level)

    def _share_long_axis(self, share, side, level):
        """When sharing main subplots, shares the long axes of their side panels,
        assuming long axis sharing is enabled for that panel."""
        # TODO: Re-calculate share settings at draw time!
        if self._panel_side:
            return
        s = side[0]
        if s not in 'lrbt':
            raise ValueError(f'Invalid side {side!r}.')
        paxs = getattr(self, '_' + s + 'panels') # calling this means, share properties on this axes with input 'share' axes
        if not all(not pax._panel_filled and pax._panel_share for pax in paxs):
            return
        axis = 'x' if s in 'tb' else 'y'
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share, level)

    def _sharex_setup(self, sharex, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharex is None or sharex is self:
            return
        if isinstance(self, ProjectionAxes) or isinstance(sharex, ProjectionAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels
        self._share_short_axis(sharex, 'l', level)
        self._share_short_axis(sharex, 'r', level)
        self._share_long_axis(sharex,  'b', level)
        self._share_long_axis(sharex,  't', level)
        # Builtin features
        if level > 0:
            self._sharex = sharex
        if level > 1:
            self._shared_x_axes.join(self, sharex)
        # "Shared" axis and tick labels
        # TODO: Why does this work?! Is only called on initialization, but
        # shouldn't this be overridden when user changes e.g. the formatter,
        # since the tick label objects themselves will also change? Maybe
        # new tick labels inherit properties from old tick labels.
        if level > 0:
            self.xaxis.label.set_visible(False)
        if level > 2:
            self.xaxis.set_major_formatter(mticker.NullFormatter())

    def _sharey_setup(self, sharey, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharey is None or sharey is self:
            return
        if isinstance(self, ProjectionAxes) or isinstance(sharey, ProjectionAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels
        self._share_short_axis(sharey, 'b', level)
        self._share_short_axis(sharey, 't', level)
        self._share_long_axis(sharey,  'l', level)
        self._share_long_axis(sharey,  'r', level)
        # Builtin features
        if level > 0:
            self._sharey = sharey
        if level > 1:
            self._shared_y_axes.join(self, sharey)
        # "Shared" axis and tick labels
        if level > 0:
            self.yaxis.label.set_visible(False)
        if level > 2:
            self.yaxis.set_major_formatter(mticker.NullFormatter())

    def _share_panels_setup(self):
        """Sets up axis sharing between main subplots and panels."""
        share = lambda paxs: (paxs and all(pax._panel_share for pax in paxs))
        # Top and bottom
        bottom = None
        if share(self._bpanels):
            bottom = self._bpanels[-1]
            for iax in (self, *self._bpanels[:-1]):
                iax._sharex_setup(bottom, 3) # parent is *bottom-most* panel
        if share(self._tpanels):
            bottom = bottom or self
            for iax in self._tpanels:
                iax._sharex_setup(bottom, 3)
        # Left and right
        left = None
        if share(self._lpanels):
            left = self._lpanels[0]
            for iax in (*self._lpanels[1:], self):
                iax._sharey_setup(left, 3) # parent is *bottom-most* panel
        if share(self._rpanels):
            left = left or self
            for iax in self._rpanels:
                iax._sharey_setup(left, 3)

    def _update_title(self, obj, **kwargs):
        """Redraws title if updating with input keyword args failed."""
        # Try to just return updated object, redraw may be necessary
        # WARNING: Making text instances invisible seems to mess up tight
        # bounding box calculations and cause other issues. Just reset text.
        keys = ('border','lw','linewidth','bordercolor','invert')
        kwextra = {key:value for key,value in kwargs.items() if key in keys}
        kwargs = {key:value for key,value in kwargs.items() if key not in keys}
        obj.update(kwargs)
        if kwextra:
            obj.set_text('')
        else:
            return obj
        # Get properties from old object
        for key in ('ha', 'va', 'color', 'transform', 'fontproperties'):
            kwextra[key] = getattr(obj, 'get_' + key)() # copy over attrs
        text = kwargs.pop('text', obj.get_text())
        x, y = kwargs.pop('position', (None,None))
        pos = obj.get_position()
        x = _notNone(kwargs.pop('x', x), pos[0])
        y = _notNone(kwargs.pop('y', y), pos[1])
        return self.text(x, y, text, **kwextra)

    def context(self, *, mode=2, rc_kw=None, **kwargs):
        """
        For internal use. Sets up temporary `~proplot.rctools.rc` settings by
        returning the result of `~proplot.rctools.rc_configurator.context`.

        Parameters
        ----------
        rc_kw : dict, optional
            A dictionary containing "rc" configuration settings that will
            be applied to this axes. Temporarily updates the
            `~proplot.rctools.rc` object. See `~proplot.rctools` for details.
        **kwargs
            Any of three options:

            * A keyword arg for `Axes.format`, `CartesianAxes.format`,
              or `ProjectionAxes.format`.
            * A global "rc" keyword arg, like ``linewidth`` or ``color``.
            * A standard "rc" keyword arg **with the dots omitted**,
              like ``landcolor`` instead of ``land.color``.

            The latter two options update the `~proplot.rctools.rc`
            object, just like `rc_kw`.

        Other parameters
        ----------------
        mode : int, optional
            The "getitem mode". This is used under-the-hood -- you shouldn't
            have to use it directly. Determines whether queries to the
            `~proplot.rctools.rc` object will ignore `rcParams <https://matplotlib.org/users/customizing.html>`__.
            This can help prevent a massive number of unnecessary lookups
            when the settings haven't been changed by the user.
            See `~proplot.rctools.rc_configurator` for details.

        Returns
        -------
        `~proplot.rctools.rc_configurator`
            The `proplot.rctools.rc` object primed for use in a "with"
            statement.
        dict
            Dictionary of keyword arguments that are not `~proplot.rctools.rc`
            properties, to be passed to the ``format`` methods.
        """
        # Figure out which kwargs are valid rc settings
        # TODO: Support for 'small', 'large', etc. font
        kw = {} # for format
        rc_kw = rc_kw or {}
        for key,value in kwargs.items():
            key_fixed = RC_NAMES_NODOTS.get(key, None)
            if key_fixed is None:
                kw[key] = value
            else:
                rc_kw[key_fixed] = value
        rc._getitem_mode = 0 # might still be non-zero if had error
        # Return "context object", which is just the configurator itself
        # primed for use in a "with" statement
        return rc.context(rc_kw, mode=mode), kw

    def format(self, *, title=None, top=None,
        figtitle=None, suptitle=None, rowlabels=None, collabels=None,
        leftlabels=None, rightlabels=None, toplabels=None, bottomlabels=None,
        llabels=None, rlabels=None, tlabels=None, blabels=None,
        **kwargs,
        ):
        """
        Called by `CartesianAxes.format` and `ProjectionAxes.format`,
        formats the axes titles, a-b-c labelling, row and column labels, and
        figure title.

        Note that the `abc`, `abcformat`, `abcloc`, and `titleloc` keyword
        arguments are actually rc configuration settings that are temporarily
        changed by the call to `~Axes.format`. They are documented here
        because it is extremely common to change them with `~Axes.format`.
        They also appear in the tables in the `~proplot.rctools` documention.

        Parameters
        ----------
        title : str, optional
            The axes title.
        abc : bool, optional
            Whether to apply "a-b-c" subplot labelling based on the
            ``number`` attribute. If ``number`` is >26, the labels will loop
            around to a, ..., z, aa, ..., zz, aaa, ..., zzz, ... God help you
            if you ever need that many labels. Defaults to ``rc['abc']``.
        abcformat : str, optional
            It is a string containing the character ``a`` or ``A``, specifying
            the format of a-b-c labels.  ``'a'`` is the default, but e.g.
            ``'a.'``, ``'a)'``, or ``'A'`` might be desirable. Defaults to
            ``rc['abc.format']``.
        abcloc, titleloc : str, optional
            They are strings indicating the location for the a-b-c label and
            main title. The following locations keys are valid. Defaults to
            ``rc['abc.loc']`` and ``rc['title.loc']``.

            ========================  ============================
            Location                  Valid keys
            ========================  ============================
            center above axes         ``'center'``, ``'c'``
            left above axes           ``'left'``, ``'l'``
            right above axes          ``'right'``, ``'r'``
            lower center inside axes  ``'lower center``', ``'lc'``
            upper center inside axes  ``'upper center'``, ``'uc'``
            upper right inside axes   ``'upper right'``, ``'ur'``
            upper left inside axes    ``'upper left'``, ``'ul'``
            lower left inside axes    ``'lower left'``, ``'ll'``
            lower right inside axes   ``'lower right'``, ``'lr'``
            ========================  ============================

        abcborder, titleborder : bool, optional
            Whether to draw a white border around titles and a-b-c labels
            positioned inside the axes. This can help them stand out on top
            of artists plotted inside the axes. Defaults to
            ``rc['abc.border']`` and ``rc['title.border']``
        ltitle, rtitle, ultitle, uctitle, urtitle, lltitle, lctitle, lrtitle : str, optional
            Axes titles in particular positions. This lets you specify multiple
            "titles" for each subplots. See the `abcloc` keyword.
        top : bool, optional
            Whether to try to put title and a-b-c label above the top subplot
            panel (if it exists), or to always put them on the main subplot.
            Defaults to ``True``.
        rowlabels, colllabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`.
        llabels, tlabels, rlabels, blabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`, `rightlabels`, `bottomlabels`.
        leftlabels, toplabels, rightlabels, bottomlabels : list of str, optional
            The subplot row and column labels. If list, length must match
            the number of subplots on the left, top, right, or bottom edges
            of the figure.
        figtitle, suptitle : str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.
            This is an improvement on matplotlib's "super" title, which just
            centers the text between figure edges.
        """
        # Figure patch (for some reason needs to be re-asserted even if
        # declared before figure is drawn)
        kw = rc.fill({'facecolor':'figure.facecolor'})
        self.figure.patch.update(kw)
        if top is not None:
            self._title_above_panel = top
        pad = rc['axes.titlepad']
        if pad is not None:
            self._set_title_offset_trans(pad)
            self._title_pad = pad

        # Super title
        # NOTE: These are actually *figure-wide* settings, but that line seems
        # to get blurred -- where we have shared axes, spanning labels, and
        # whatnot. May result in redundant assignments if formatting more than
        # one axes, but operations are fast so some redundancy is nbd.
        fig = self.figure
        suptitle = _notNone(figtitle, suptitle, None, names=('figtitle','suptitle'))
        kw = rc.fill({
            'fontsize':   'suptitle.fontsize',
            'weight':     'suptitle.weight',
            'color':      'suptitle.color',
            'fontfamily': 'font.family'
            })
        if suptitle or kw:
            fig._update_suptitle(suptitle, **kw)
        # Labels
        llabels = _notNone(rowlabels, leftlabels, llabels, None, names=('rowlabels','leftlabels','llabels'))
        tlabels = _notNone(collabels, toplabels, tlabels, None, names=('collabels','toplabels','tlabels'))
        rlabels = _notNone(rightlabels, rlabels, None, names=('rightlabels','rlabels'))
        blabels = _notNone(bottomlabels, blabels, None, names=('bottomlabels','blabels'))
        for side,labels in zip(
            ('left', 'right', 'top', 'bottom'),
            (llabels, rlabels, tlabels, blabels),
            ):
            kw = rc.fill({
                'fontsize':   side + 'label.fontsize',
                'weight':     side + 'label.weight',
                'color':      side + 'label.color',
                'fontfamily': 'font.family'
                })
            if labels or kw:
                fig._update_suplabels(self, side, labels, **kw)

        # A-b-c labels
        titles_dict = self._titles_dict
        if not self._panel_side:
            # Location and text
            abcformat = rc['abc.format'] # changed, or running format for first time?
            if abcformat and self.number is not None:
                if 'a' not in abcformat and 'A' not in abcformat:
                    raise ValueError(f'Invalid abcformat {abcformat!r}. Must include letter "a" or "A".')
                abcedges = abcformat.split('a' if 'a' in abcformat else 'A')
                text = abcedges[0] + _abc(self.number-1) + abcedges[-1]
                if 'A' in abcformat:
                    text = text.upper()
                self._abc_text = text
            # Apply new settings
            # Also if a-b-c label was moved, remove previous one and update
            # text on new one, in case self._abc_text has not changed.
            loc, obj, kw = self._get_title_props(abc=True)
            iloc = self._abc_loc
            obj = self._update_title(obj, **kw)
            titles_dict[loc] = obj
            if iloc is not None and loc!=iloc:
                self.abc.set_text('')
                obj.set_text(self._abc_text)
            self.abc = obj
            self._abc_loc = loc
            # Toggle visibility
            # NOTE: If abc is a matplotlib 'title' attribute, making it
            # invisible messes stuff up. Just set text to empty.
            abc = rc['abc']
            if abc is not None:
                obj.set_text(self._abc_text if bool(abc) else '')

        # Titles
        # Tricky because we have to reconcile two workflows:
        # 1. title='name' and titleloc='position'
        # 2. ltitle='name', rtitle='name', etc., arbitrarily many titles
        # First update existing titles
        # NOTE: _update_title should never return new objects unless called
        # with *inner* titles... *outer* titles will just refresh, so we
        # don't need to re-assign the attributes or anything.
        loc, obj, kw = self._get_title_props()
        if kw:
            for iloc,iobj in titles_dict.items():
                if iloc is self._abc_loc:
                    continue
                titles_dict[iloc] = self._update_title(iobj, **kw)
        # Workflow 2, want this to come first so workflow 1 gets priority
        for ikey,ititle in kwargs.items():
            if not ikey[-5:] == 'title':
                raise ValueError(f'format() got an unexpected keyword argument {ikey!r}.')
            iloc, iobj, ikw = self._get_title_props(loc=ikey[:-5])
            if ititle is not None:
                ikw['text'] = ititle
            if ikw:
                titles_dict[iloc] = self._update_title(iobj, **ikw)
        # Workflow 1, make sure that if user calls ax.format(title='Title')
        # *then* ax.format(titleloc='left') it copies over the text.
        iloc = self._title_loc
        if iloc is not None and loc != iloc:
            iobj = titles_dict[iloc]
            if title is None:
                title = iobj.get_text()
            iobj.set_text('')
        self._title_loc = loc # assigns default loc on first run, or changed loc
        if title is not None:
            kw['text'] = title
        if kw:
            titles_dict[loc] = self._update_title(obj, **kw)

    def area(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.fill_between`, which is wrapped by
        `~proplot.wrappers.fill_between_wrapper`."""
        # NOTE: *Cannot* assign area = axes.Axes.fill_between because the
        # wrapper won't be applied and for some reason it messes up
        # automodsumm, which tries to put the matplotlib docstring on website
        return self.fill_between(*args, **kwargs)

    def areax(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.fill_betweenx`, which is wrapped by
        `~proplot.wrappers.fill_betweenx_wrapper`."""
        return self.fill_betweenx(*args, **kwargs)

    def boxes(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.boxplot`, which is wrapped by
        `~proplot.wrappers.boxplot_wrapper`."""
        return self.boxplot(*args, **kwargs)

    def cmapline(self, *args, values=None,
        cmap=None, norm=None,
        interp=0, **kwargs):
        """
        Invoked by `~proplot.wrappers.plot_wrapper` when you pass the `cmap`
        keyword argument to `~matplotlib.axes.Axes.plot`. Draws a "colormap line",
        i.e. a line whose color changes as a function of some parametric coordinate
        `values`. This is actually a collection of lines, added as a
        `~matplotlib.collections.LineCollection` instance. See `this matplotlib example
        <https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html>`__.

        Parameters
        ----------
        *args : (y,) or (x,y)
            The coordinates. If `x` is not provided, it is inferred from `y`.
        cmap : colormap spec, optional
            The colormap specifier, passed to `~proplot.styletools.Colormap`.
        values : list of float
            The parametric values used to map points on the line to colors
            in the colormap.
        norm : normalizer spec, optional
            The normalizer, passed to `~proplot.styletools.Norm`.
        interp : int, optional
            Number of values between each line joint and each *halfway* point
            between line joints to which you want to interpolate.
        """
        # First error check
        # WARNING: So far this only works for 1D *x* and *y* coordinates. Cannot
        # draw multiple colormap lines at once, unlike `~matplotlib.axes.Axes.plot`.
        if values is None:
            raise ValueError('Requires a "values" keyword arg.')
        if len(args) not in (1,2):
            raise ValueError(f'Requires 1-2 arguments, got {len(args)}.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(y.shape[-1]) if len(args) == 1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim != 1 or y.ndim != 1 or values.ndim != 1:
            raise ValueError(f'x ({x.ndim}-d), y ({y.ndim}-d), and values ({values.ndim}-d) must be 1-dimensional.')
        if len(x) != len(y) or len(x) != len(values) or len(y) != len(values):
            raise ValueError(f'{len(x)} xs, {len(y)} ys, but {len(values)} colormap values.')

        # Interpolate values to allow for smooth gradations between values
        # (bins=False) or color switchover halfway between points (bins=True)
        # Then optionally interpolate the corresponding colormap values
        if interp > 0:
            xorig, yorig, vorig = x, y, values
            x, y, values = [], [], []
            for j in range(xorig.shape[0]-1):
                idx = (slice(None, -1) if j+1 < xorig.shape[0]-1 else slice(None))
                x.extend(np.linspace(xorig[j], xorig[j+1], interp + 2)[idx].flat)
                y.extend(np.linspace(yorig[j], yorig[j+1], interp + 2)[idx].flat)
                values.extend(np.linspace(vorig[j], vorig[j+1], interp + 2)[idx].flat)
            x, y, values = np.array(x), np.array(y), np.array(values)
        coords = []
        levels = utils.edges(values)
        for j in range(y.shape[0]):
            # Get x/y coordinates and values for points to the 'left' and
            # 'right' of each joint
            if j == 0:
                xleft, yleft = [], []
            else:
                xleft = [(x[j-1] + x[j])/2, x[j]]
                yleft = [(y[j-1] + y[j])/2, y[j]]
            if j+1 == y.shape[0]:
                xright, yright = [], []
            else:
                xleft  = xleft[:-1] # prevent repetition when joined with right
                yleft  = yleft[:-1]
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

    def colorbar(self, *args, loc=None, pad=None,
        length=None, width=None, space=None, frame=None, frameon=None,
        alpha=None, linewidth=None, edgecolor=None, facecolor=None,
        **kwargs):
        """
        Adds colorbar as an *inset* or along the outside edge of the axes.
        See `~proplot.wrappers.colorbar_wrapper` for details.

        Parameters
        ----------
        loc : str, optional
            The colorbar location. Defaults to ``rc['colorbar.loc']``. The
            following location keys are valid.

            ==================  ==========================================================
            Location            Valid keys
            ==================  ==========================================================
            outer left          ``'l'``, ``'left'``
            outer right         ``'r'``, ``'right'``
            outer bottom        ``'b'``, ``'bottom'``
            outer top           ``'t'``, ``'top'``
            default inset       ``0``, ``'i'``, ``'inset'``
            upper right inset   ``1``, ``'upper right'``, ``'ur'``
            upper left inset    ``2``, ``'upper left'``, ``'ul'``
            lower left inset    ``3``, ``'lower left'``, ``'ll'``
            lower right inset   ``4``, ``'lower right'``, ``'lr'``
            ==================  ==========================================================

        pad : float or str, optional
            The space between the axes edge and the colorbar. Ignored for
            outer colorbars. If float, units are inches. If string, units
            are interpreted by `~proplot.utils.units`. Defaults to
            ``rc['colorbar.pad']``.
        length : float or str, optional
            The colorbar length. For outer colorbars, units are relative to
            the axes width or height. For inset colorbars, if float, units are
            inches; if string, units are interpreted by `~proplot.utils.units`.
            Defaults to ``rc['colorbar.length']`` for outer colorbars,
            ``rc['colorbar.lengthinset']`` for inset colorbars.
        width : float or str, optional
            The colorbar width. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. Defaults to
            ``rc['colorbar.width']`` for outer colorbars,
            ``rc['colorbar.widthinset']`` for inset colorbars.
        space : float or str, optional
            The space between the colorbar and the main axes for outer
            colorbars. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. By default, this
            is adjusted automatically in the "tight layout" calculation, or is
            ``rc['subplots.panelspace']`` if "tight layout" is turned off.
        frame, frameon : bool, optional
            Whether to draw a frame around inset colorbars, just like
            `~matplotlib.axes.Axes.legend`.
            Defaults to ``rc['colorbar.frameon']``.
        alpha, linewidth, edgecolor, facecolor : optional
            Transparency, edge width, edge color, and face color for the frame
            around the inset colorbar. Defaults to
            ``rc['colorbar.framealpha']``, ``rc['axes.linewidth']``,
            ``rc['axes.edgecolor']``, and ``rc['axes.facecolor']``,
            respectively.
        **kwargs
            Passed to `~proplot.wrappers.colorbar_wrapper`.
        """
        # TODO: add option to pad inset away from axes edge!
        kwargs.update({'edgecolor':edgecolor, 'linewidth':linewidth})
        loc = _notNone(loc, rc['colorbar.loc'])
        loc = self._loc_translate(loc)
        if loc == 'best': # a white lie
            loc = 'lower right'
        if not isinstance(loc, str): # e.g. 2-tuple or ndarray
            raise ValueError(f'Invalid colorbar location {loc!r}.')

        # Generate panel
        if loc in ('left','right','top','bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            return ax.colorbar(loc='_fill', *args, **kwargs)

        # Filled colorbar
        if loc == '_fill':
            # Hide content and resize panel
            # NOTE: Do not run self.clear in case we want title above this
            for s in self.spines.values():
                s.set_visible(False)
            self.xaxis.set_visible(False)
            self.yaxis.set_visible(False)
            self.patch.set_alpha(0)
            self._panel_filled = True
            # Draw colorbar with arbitrary length relative to full length of panel
            side = self._panel_side
            length = _notNone(length, rc['colorbar.length'])
            subplotspec = self.get_subplotspec()
            if length <= 0 or length > 1:
                raise ValueError(f'Panel colorbar length must satisfy 0 < length <= 1, got length={length!r}.')
            if side in ('bottom','top'):
                gridspec = mgridspec.GridSpecFromSubplotSpec(
                        nrows=1, ncols=3, wspace=0,
                        subplot_spec=subplotspec,
                        width_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subplotspec = gridspec[1]
            else:
                gridspec = mgridspec.GridSpecFromSubplotSpec(
                        nrows=3, ncols=1, hspace=0,
                        subplot_spec=subplotspec,
                        height_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subplotspec = gridspec[1]
            with self.figure._unlock():
                ax = self.figure.add_subplot(subplotspec, projection=None)
            if ax is self:
                raise ValueError(f'Uh oh.')
            self.add_child_axes(ax)
            # Location
            if side in ('bottom','top'):
                outside, inside = 'bottom', 'top'
                if side == 'top':
                    outside, inside = inside, outside
                ticklocation = outside
                orientation  = 'horizontal'
            else:
                outside, inside = 'left', 'right'
                if side == 'right':
                    outside, inside = inside, outside
                ticklocation = outside
                orientation  = 'vertical'
            # Keyword args and add as child axes
            orient = kwargs.get('orientation', None)
            if orient is not None and orient != orientation:
                warnings.warn(f'Overriding input orientation={orient!r}.')
            ticklocation = kwargs.pop('tickloc', None) or ticklocation
            ticklocation = kwargs.pop('ticklocation', None) or ticklocation
            kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})

        # Inset colorbar
        else:
            # Default props
            cbwidth, cblength = width, length
            width, height = self.get_size_inches()
            extend = units(_notNone(kwargs.get('extendsize',None), rc['colorbar.extendinset']))
            cbwidth = units(_notNone(cbwidth, rc['colorbar.widthinset']))/height
            cblength = units(_notNone(cblength, rc['colorbar.lengthinset']))/width
            pad = units(_notNone(pad, rc['colorbar.axespad']))
            xpad, ypad = pad/width, pad/height
            # Get location in axes-relative coordinates
            # Bounds are x0, y0, width, height in axes-relative coordinate to start
            if kwargs.get('label', ''):
                xspace = 2.4*rc['font.size']/72 + rc['xtick.major.size']/72
            else:
                xspace = 1.2*rc['font.size']/72 + rc['xtick.major.size']/72
            xspace /= height # space for labels
            if loc == 'upper right':
                bounds = (1 - xpad - cblength, 1 - ypad - cbwidth)
                fbounds = (1 - 2*xpad - cblength, 1 - 2*ypad - cbwidth - xspace)
            elif loc == 'upper left':
                bounds = (xpad, 1 - ypad - cbwidth)
                fbounds = (0, 1 - 2*ypad - cbwidth - xspace)
            elif loc == 'lower left':
                bounds = (xpad, ypad + xspace)
                fbounds = (0, 0)
            elif loc == 'lower right':
                bounds = (1 - xpad - cblength, ypad + xspace)
                fbounds = (1 - 2*xpad - cblength, 0)
            else:
                raise ValueError(f'Invalid colorbar location {loc!r}.')
            bounds = (bounds[0], bounds[1], cblength, cbwidth)
            fbounds = (fbounds[0], fbounds[1], 2*xpad + cblength, 2*ypad + cbwidth + xspace)
            # Make frame
            # NOTE: We do not allow shadow effects or fancy edges effect.
            # Also keep zorder same as with legend.
            frameon = _notNone(frame, frameon, rc['colorbar.frameon'], names=('frame','frameon'))
            if frameon:
                # Make patch object
                xmin, ymin, width, height = fbounds
                patch = mpatches.Rectangle((xmin,ymin), width, height,
                        snap=True, zorder=4.5, transform=self.transAxes)
                # Update patch props
                alpha = _notNone(alpha, rc['colorbar.framealpha'])
                linewidth = _notNone(linewidth, rc['axes.linewidth'])
                edgecolor = _notNone(edgecolor, rc['axes.edgecolor'])
                facecolor = _notNone(facecolor, rc['axes.facecolor'])
                patch.update({'alpha':alpha, 'linewidth':linewidth,
                              'edgecolor':edgecolor, 'facecolor':facecolor})
                self.add_artist(patch)
            # Make axes
            locator = self._make_inset_locator(bounds, self.transAxes)
            bbox = locator(None, None)
            ax = maxes.Axes(self.figure, bbox.bounds, zorder=5)
            ax.set_axes_locator(locator)
            self.add_child_axes(ax)
            # Default keyword args
            orient = kwargs.pop('orientation', None)
            if orient is not None and orient != 'horizontal':
                warnings.warn(f'Orientation for inset colorbars must be horizontal, ignoring orient={orient!r}.')
            ticklocation = kwargs.pop('tickloc', None)
            ticklocation = kwargs.pop('ticklocation', None) or ticklocation
            if ticklocation is not None and ticklocation != 'bottom':
                warnings.warn(f'Inset colorbars can only have ticks on the bottom.')
            kwargs.update({'orientation':'horizontal', 'ticklocation':'bottom'})
            kwargs.setdefault('maxn', 5)
            kwargs.setdefault('extendsize', extend)

        # Generate colorbar
        return wrappers.colorbar_wrapper(ax, *args, **kwargs)

    def legend(self, *args, loc=None, width=None, space=None, **kwargs):
        """
        Adds an *inset* legend or *outer* legend along the edge of the axes.
        See `~matplotlib.axes.Axes.legend` and
        `~proplot.wrappers.legend_wrapper` for details.

        Parameters
        ----------
        loc : int or str, optional
            The legend location or panel location. The following location keys
            are valid. Note that if a panel does not exist, it will be
            generated on-the-fly.

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            left panel          ``'l'``, ``'left'``
            right panel         ``'r'``, ``'right'``
            bottom panel        ``'b'``, ``'bottom'``
            top panel           ``'t'``, ``'top'``
            "best" inset        ``0``, ``'best'``, ``'inset'``, ``'i'``
            upper right inset   ``1``, ``'upper right'``, ``'ur'``
            upper left inset    ``2``, ``'upper left'``, ``'ul'``
            lower left inset    ``3``, ``'lower left'``, ``'ll'``
            lower right inset   ``4``, ``'lower right'``, ``'lr'``
            center left inset   ``5``, ``'center left'``, ``'cl'``
            center right inset  ``6``, ``'center right'``, ``'cr'``
            lower center inset  ``7``, ``'lower center'``, ``'lc'``
            upper center inset  ``8``, ``'upper center'``, ``'uc'``
            center inset        ``9``, ``'center'``, ``'c'``
            ==================  =======================================

        width : float or str, optional
            The space allocated for the outer legends. This does nothing
            if "tight layout" is turned on. If float, units are inches. If
            string, units are interpreted by `~proplot.utils.units`.
        space : float or str, optional
            The space between the axes and the legend for outer legends.
            If float, units are inches. If string, units are interpreted by
            `~proplot.utils.units`.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~matplotlib.axes.Axes.legend`.
        """
        loc = self._loc_translate(loc, width=width, space=space)
        if isinstance(loc, np.ndarray):
            loc = loc.tolist()

        # Generate panel
        if loc in ('left','right','top','bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            return ax.legend(*args, loc='_fill', **kwargs)

        # Fill
        if loc == '_fill':
            # Hide content
            for s in self.spines.values():
                s.set_visible(False)
            self.xaxis.set_visible(False)
            self.yaxis.set_visible(False)
            self.patch.set_alpha(0)
            self._panel_filled = True
            # Try to make handles and stuff flush against the axes edge
            kwargs.setdefault('borderaxespad', 0)
            frameon = _notNone(kwargs.get('frame', None), kwargs.get('frameon', None), rc['legend.frameon'])
            if not frameon:
                kwargs.setdefault('borderpad', 0)
            # Apply legend location
            side = self._panel_side
            if side == 'bottom':
                loc = 'upper center'
            elif side == 'right':
                loc = 'center left'
            elif side == 'left':
                loc = 'center right'
            elif side == 'top':
                loc = 'lower center'
            else:
                raise ValueError(f'Invalid panel side {side!r}.')

        # Draw legend
        return wrappers.legend_wrapper(self, *args, loc=loc, **kwargs)

    def draw(self, renderer=None, *args, **kwargs):
        """Adds post-processing steps before axes is drawn."""
        self._reassign_title()
        super().draw(renderer, *args, **kwargs)

    def get_size_inches(self):
        """Returns the width and the height of the axes in inches."""
        width, height = self.figure.get_size_inches()
        width = width*abs(self.get_position().width)
        height = height*abs(self.get_position().height)
        return width, height

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Adds post-processing steps before tight bounding box is
        calculated, and stores the bounding box as an attribute."""
        self._reassign_title()
        bbox = super().get_tightbbox(renderer, *args, **kwargs)
        self._tight_bbox = bbox
        return bbox

    def heatmap(self, *args, **kwargs):
        """Calls `~matplotlib.axes.Axes.pcolormesh` and applies default formatting
        that is suitable for heatmaps: no gridlines, no minor ticks, and major
        ticks at the center of each grid box."""
        obj = self.pcolormesh(*args, **kwargs)
        xlocator, ylocator = None, None
        if hasattr(obj, '_coordinates'): # be careful in case private API changes! but this is only way to infer coordinates
            coords = obj._coordinates
            coords = (coords[1:,...] + coords[:-1,...])/2
            coords = (coords[:,1:,:] + coords[:,:-1,:])/2
            xlocator, ylocator = coords[0,:,0], coords[:,0,1]
        self.format(
            xgrid=False, ygrid=False, xtickminor=False, ytickminor=False,
            xlocator=xlocator, ylocator=ylocator,
            )
        return obj

    def inset_axes(self, bounds, *, transform=None, zorder=5,
        zoom=True, zoom_kw=None, **kwargs):
        """
        Like the builtin `~matplotlib.axes.Axes.inset_axes` method, but
        draws an inset `CartesianAxes` axes and adds some options.

        Parameters
        ----------
        bounds : list of float
            The bounds for the inset axes, listed as ``(x, y, width, height)``.
        transform : {'data', 'axes', 'figure'} or `~matplotlib.transforms.Transform`, optional
            The transform used to interpret `bounds`. Can be a
            `~matplotlib.transforms.Transform` object or a string representing the
            `~matplotlib.axes.Axes.transData`, `~matplotlib.axes.Axes.transAxes`,
            or `~matplotlib.figure.Figure.transFigure` transforms. Defaults to
            ``'axes'``, i.e. `bounds` is in axes-relative coordinates.
        zorder : float, optional
            The zorder of the axes, should be greater than the zorder of
            elements in the parent axes. Defaults to ``5``.
        zoom : bool, optional
            Whether to draw lines indicating the inset zoom using
            `~Axes.indicate_inset_zoom`. The lines will automatically
            adjust whenever the parent axes or inset axes limits are changed.
            Defaults to ``True``.
        zoom_kw : dict, optional
            Passed to `~Axes.indicate_inset_zoom`.

        Other parameters
        ----------------
        **kwargs
            Passed to `CartesianAxes`.
        """
        # Carbon copy with my custom axes
        if not transform:
            transform = self.transAxes
        else:
            transform = wrappers._get_transform(self, transform)
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        locator = self._make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = CartesianAxes(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)
        # The following locator lets the axes move if we used data coordinates,
        # is called by ax.apply_aspect()
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        ax._inset_zoom = zoom
        ax._inset_parent = self
        # Zoom indicator (NOTE: Requires version >=3.0)
        if zoom:
            zoom_kw = zoom_kw or {}
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def indicate_inset_zoom(self, alpha=None,
        lw=None, linewidth=None,
        color=None, edgecolor=None, **kwargs):
        """
        Called automatically when using `~Axes.inset` with ``zoom=True``.
        Like `~matplotlib.axes.Axes.indicate_inset_zoom`, but *refreshes* the
        lines at draw-time.

        This method is called from the *inset* axes, not the parent axes.

        Parameters
        ----------
        alpha : float, optional
            The transparency of the zoom box fill.
        lw, linewidth : float, optional
            The width of the zoom lines and box outline in points.
        color, edgecolor : color-spec, optional
            The color of the zoom lines and box outline.
        **kwargs
            Passed to `~matplotlib.axes.Axes.indicate_inset`.
        """
        # Should be called from the inset axes
        parent = self._inset_parent
        alpha = alpha or 1.0
        linewidth = _notNone(lw, linewidth, rc['axes.linewidth'], names=('lw', 'linewidth'))
        edgecolor = _notNone(color, edgecolor, rc['axes.edgecolor'], names=('color', 'edgecolor'))
        if not parent:
            raise ValueError(f'{self} is not an inset axes.')
        xlim, ylim = self.get_xlim(), self.get_ylim()
        rect = (xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0])
        # Call indicate_inset
        rectpatch, connects = parent.indicate_inset(rect, self,
            linewidth=linewidth, edgecolor=edgecolor, alpha=alpha, **kwargs)
        # Adopt properties from old one
        if self._inset_zoom_data:
            rectpatch_old, connects_old = self._inset_zoom_data
            rectpatch.update_from(rectpatch_old)
            rectpatch_old.set_visible(False)
            for line,line_old in zip(connects,connects_old):
                visible = line.get_visible()
                line.update_from(line_old)
                line.set_linewidth(line_old.get_linewidth())
                line.set_visible(visible)
                line_old.set_visible(False)
        # Format zoom data
        else:
            for line in connects:
                line.set_linewidth(linewidth)
                line.set_color(edgecolor)
                line.set_alpha(alpha)
        self._inset_zoom_data = (rectpatch, connects)
        return (rectpatch, connects)

    def panel_axes(self, side, **kwargs):
        """
        Returns a panel drawn along the edge of an axes.

        Parameters
        ----------
        ax : `~proplot.axes.Axes`
            The axes for which we are drawing a panel.
        width : float or str or list thereof, optional
            The panel width. If float, units are inches. If string, units are
            interpreted by `~proplot.utils.units`.
        space : float or str or list thereof, optional
            Empty space between the main subplot and the panel. If float,
            units are inches. If string, units are interpreted by
            `~proplot.utils.units`.
        share : bool, optional
            Whether to enable axis sharing between the *x* and *y* axes of the
            main subplot and the panel long axes for each panel in the stack.
            Sharing between the panel short axis and other panel short axes
            is determined by figure-wide `sharex` and `sharey` settings.

        Returns
        -------
        `~proplot.axes.Axes`
            The panel axes.
        """
        return self.figure._add_axes_panel(self, side, **kwargs)

    def violins(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.violinplot`, which is wrapped by
        `~proplot.wrappers.violinplot_wrapper`."""
        return self.violinplot(*args, **kwargs)

    panel = panel_axes
    """Alias for `~Axes.panel_axes`."""
    inset = inset_axes
    """Alias for `~Axes.inset_axes`."""

    @property
    def number(self):
        """The axes number, controls a-b-c label order and order of
        appearence in the `~proplot.subplots.axes_grid` returned by
        `~proplot.subplots.subplots`."""
        return self._number

    def _iter_panels(self, sides='lrbt'):
        """Iterates over axes and child panel axes."""
        axs = [self] if self.get_visible() else []
        if not ({*sides} <= {*'lrbt'}):
            raise ValueError(f'Invalid sides {sides!r}.')
        for s in sides:
            for ax in getattr(self, '_' + s + 'panels'):
                if not ax or not ax.get_visible():
                    continue
                axs.append(ax)
        return axs

#-----------------------------------------------------------------------------#
# Axes subclasses
#-----------------------------------------------------------------------------#
def _rcloc_to_stringloc(x, string): # figures out string location
    """Gets *location string* from the *boolean* rc settings, for a given
    string prefix like ``'axes.spines'`` or ``'xtick'``."""
    if x == 'x':
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
    elif x == 'y':
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
        raise ValueError(f'"x" must equal "x" or "y".')

class CartesianAxes(Axes):
    """
    Axes subclass for ordinary Cartesian axes. Adds several new methods and
    overrides existing ones.

    See also
    --------
    `~proplot.subplots.subplots`, `Axes`
    """
    name = 'cartesian'
    """The registered projection name."""
    def __init__(self, *args, **kwargs):
        # Impose default formatter
        super().__init__(*args, **kwargs)
        formatter = axistools.Formatter('default')
        self.xaxis.set_major_formatter(formatter)
        self.yaxis.set_major_formatter(formatter)
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True
        # Custom attributes
        self._datex_rotated = False # whether to automatically apply rotation to datetime labels during post processing
        self._dualy_scale = None # for scaling units on opposite side of ax
        self._dualx_scale = None

    def __getattribute__(self, attr, *args):
        """Applies the `~proplot.wrappers.cmap_wrapper`,
        `~proplot.wrappers.cycle_wrapper`, `~proplot.wrappers.add_errorbars`,
        `~proplot.wrappers.enforce_centers`, `~proplot.wrappers.enforce_edges`,
        `~proplot.wrappers.plot_wrapper`, `~proplot.wrappers.scatter_wrapper`,
        `~proplot.wrappers.bar_wrapper`, `~proplot.wrappers.barh_wrapper`,
        `~proplot.wrappers.boxplot_wrapper`, `~proplot.wrappers.violinplot_wrapper`,
        `~proplot.wrappers.fill_between_wrapper`, and `~proplot.wrappers.fill_betweenx_wrapper`
        wrappers."""
        obj = super().__getattribute__(attr, *args)
        if callable(obj):
            # Step 3) Color usage wrappers
            if attr in wrappers.CMAP_METHODS: # must come first!
                obj = wrappers._cmap_wrapper(self, obj)
            elif attr in wrappers.CYCLE_METHODS:
                obj = wrappers._cycle_wrapper(self, obj)
            # Step 2) Utilities
            if attr in wrappers.CENTERS_METHODS:
                obj = wrappers._enforce_centers(self, obj)
            elif attr in wrappers.EDGES_METHODS:
                obj = wrappers._enforce_edges(self, obj)
            elif attr in wrappers.ERRORBAR_METHODS:
                obj = wrappers._add_errorbars(self, obj)
            # Step 1) Parse input
            if attr in wrappers.D2_METHODS:
                obj = wrappers._autoformat_2d(self, obj)
            elif attr in wrappers.D1_METHODS:
                obj = wrappers._autoformat_1d(self, obj)
            # Step 0) Special wrappers
            if attr == 'plot':
                obj = wrappers._plot_wrapper(self, obj)
            elif attr == 'scatter':
                obj = wrappers._scatter_wrapper(self, obj)
            elif attr == 'boxplot':
                obj = wrappers._boxplot_wrapper(self, obj)
            elif attr == 'violinplot':
                obj = wrappers._violinplot_wrapper(self, obj)
            elif attr == 'bar':
                obj = wrappers._bar_wrapper(self, obj)
            elif attr == 'barh': # skips cycle wrapper and calls bar method
                obj = wrappers._barh_wrapper(self, obj)
            elif attr == 'hist': # skips cycle wrapper and calls bar method
                obj = wrappers._hist_wrapper(self, obj)
            elif attr == 'fill_between':
                obj = wrappers._fill_between_wrapper(self, obj)
            elif attr == 'fill_betweenx':
                obj = wrappers._fill_betweenx_wrapper(self, obj)
        return obj

    def _altx_overrides(self):
        """Applies alternate *x* axis overrides."""
        # Unlike matplotlib API, we strong arm user into certain twin axes
        # settings... doesn't really make sense to have twin axes without this
        if self._altx_child is not None: # altx was called on this axes
            self.spines['top'].set_visible(False)
            self.spines['bottom'].set_visible(True)
            self.xaxis.tick_bottom()
            self.xaxis.set_label_position('bottom')
        elif self._altx_parent is not None: # this axes is the result of altx
            self.spines['bottom'].set_visible(False)
            self.spines['top'].set_visible(True)
            self.spines['left'].set_visible(False)
            self.spines['right'].set_visible(False)
            self.xaxis.tick_top()
            self.xaxis.set_label_position('top')
            self.yaxis.set_visible(False)
            self.patch.set_visible(False)

    def _alty_overrides(self):
        """Applies alternate *y* axis overrides."""
        if self._alty_child is not None:
            self.spines['right'].set_visible(False)
            self.spines['left'].set_visible(True)
            self.yaxis.tick_left()
            self.yaxis.set_label_position('left')
        elif self._alty_parent is not None:
            self.spines['left'].set_visible(False)
            self.spines['right'].set_visible(True)
            self.spines['top'].set_visible(False)
            self.spines['bottom'].set_visible(False)
            self.yaxis.tick_right()
            self.yaxis.set_label_position('right')
            self.xaxis.set_visible(False)
            self.patch.set_visible(False)

    def _datex_rotate(self):
        """Applies default rotation to datetime axis coordinates."""
        # NOTE: Rotation is done *before* horizontal/vertical alignment,
        # cannot change alignment with set_tick_params. Must apply to text
        # objects. fig.autofmt_date calls subplots_adjust, so cannot use it.
        if (not isinstance(self.xaxis.converter, mdates.DateConverter) or
            self._datex_rotated): # user applied custom rotation
            return
        rotation = rc['axes.formatter.timerotation']
        kw = {'rotation':rotation}
        if rotation not in (0,90,-90):
            kw['ha'] = ('right' if rotation > 0 else 'left')
        for label in self.xaxis.get_ticklabels():
            label.update(kw)
        self._datex_rotated = True # do not need to apply more than once

    def _dualx_lock(self):
        """Locks child "dual" x-axis limits to the parent."""
        # Why did I copy and paste the dualx/dualy code you ask? Copy
        # pasting is bad, but so are a bunch of ugly getattr(attr)() calls
        scale = self._dualx_scale
        if scale is None:
            return
        child = self._altx_child
        xlim = self.get_xlim()
        xscale = self.get_xscale()
        if any(np.array(xlim) <= 0) and re.match('^(log|inverse)', xscale):
            raise ValueError(f'Axis limits go negative, but "alternate units" axis uses {xscale!r} scale.')
        # Transform axis limits with axis scale transform
        transform = child.xaxis._scale.get_transform()
        nlim = transform.inverted().transform(np.array(xlim))
        # If the transform flipped the limits, when we set axis limits, it
        # will get flipped again! So reverse the flip
        if np.sign(np.diff(xlim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]
        child.set_xlim(scale[0] + scale[1]*nlim)

    def _dualy_lock(self):
        """Locks child "dual" y-axis limits to the parent."""
        scale = self._dualy_scale
        if scale is None:
            return
        child = self._alty_child
        ylim = self.get_ylim()
        yscale = self.get_yscale()
        if any(np.array(ylim) <= 0) and re.match('^(log|inverse)', yscale):
            raise ValueError(f'Axis limits go negative, but "alternate units" axis uses {yscale!r} scale.')
        transform = child.yaxis._scale.get_transform()
        nlim = transform.inverted().transform(np.array(ylim))
        if np.sign(np.diff(ylim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]
        child.set_ylim(scale[0] + scale[1]*nlim)

    def format(self, *,
        aspect=None,
        xloc=None, yloc=None,
        xspineloc=None, yspineloc=None,
        xtickloc=None, ytickloc=None, fixticks=False,
        xlabelloc=None, ylabelloc=None,
        xticklabelloc=None, yticklabelloc=None,
        xtickdir=None, ytickdir=None,
        xgrid=None, ygrid=None,
        xgridminor=None, ygridminor=None,
        xtickminor=True, ytickminor=True,
        xticklabeldir=None, yticklabeldir=None,
        xtickrange=None, ytickrange=None,
        xreverse=False, yreverse=False,
        xlabel=None, ylabel=None,
        xlim=None, ylim=None,
        xscale=None, yscale=None,
        xrotation=None, yrotation=None,
        xformatter=None, yformatter=None, xticklabels=None, yticklabels=None,
        xticks=None, xminorticks=None, xlocator=None, xminorlocator=None,
        yticks=None, yminorticks=None, ylocator=None, yminorlocator=None,
        xbounds=None, ybounds=None,
        xmargin=None, ymargin=None,
        xcolor=None, ycolor=None,
        xticklen=None, yticklen=None,
        xlinewidth=None, ylinewidth=None,
        xlabel_kw=None, ylabel_kw=None,
        xscale_kw=None, yscale_kw=None,
        xlocator_kw=None, ylocator_kw=None,
        xformatter_kw=None, yformatter_kw=None,
        xminorlocator_kw=None, yminorlocator_kw=None,
        patch_kw=None,
        **kwargs):
        """
        Calls `Axes.format` and `Axes.context`, formats the
        *x* and *y* axis labels, tick locations, tick labels,
        axis scales, spine settings, and more.

        Parameters
        ----------
        aspect : {'auto', 'equal'}, optional
            The aspect ratio mode. If ``'auto'``, the aspect ratio is
            determined from the *x* and *y* axis limits, and ProPlot adjusts
            the subplot layout to remove excessive whitespace.
        xlabel, ylabel : str, optional
            The *x* and *y* axis labels. Applied with
            `~matplotlib.axes.Axes.set_xlabel`
            and `~matplotlib.axes.Axes.set_ylabel`.
        xlabel_kw, ylabel_kw : dict-like, optional
            The *x* and *y* axis label settings. Applied with the
            `~matplotlib.artist.Artist.update` method on the
            `~matplotlib.text.Text` instance. Options include ``'color'``,
            ``'size'``, and ``'weight'``.
        xlim, ylim : (float or None, float or None), optional
            The *x* and *y* axis data limits. Applied with
            `~matplotlib.axes.Axes.set_xlim` and `~matplotlib.axes.Axes.set_ylim`.
        xreverse, yreverse : bool, optional
            Puts the *x* and *y* axis in reverse, i.e. going from positive to
            negative numbers.
        xscale, yscale : axis scale spec, optional
            The *x* and *y* axis scales. Passed to the `~proplot.axistools.Scale`
            constructor. For example, ``xscale='log'`` applies logarithmic
            scaling, and ``xscale=('cutoff', 0.5, 2)`` applies scaling according
            to the class generated by `plot.CutoffScaleFactory('cutoff', 0.5, 2)`.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.axistools.Scale`.
        xspineloc, yspineloc : {'both', 'bottom', 'top', 'left', 'right', 'neither', 'center', 'zero'}, optional
            The *x* and *y* axis spine locations.
        xloc, yloc : optional
            Aliases for `xspineloc`, `yspineloc`.
        xtickloc, ytickloc : {'both', 'bottom', 'top', 'left', 'right', 'neither'}, optional
            Which *x* and *y* axis spines should have major and minor tick
            marks.
        xtickminor, ytickminor : bool, optional
            Whether to draw minor ticks on the *x* and *y* axes.
        xtickdir, ytickdir : {'out', 'in', 'inout'}
            Direction that major and minor tick marks point for the *x* and
            *y* axis.
        xgrid, ygrid : bool, optional
            Whether to draw major gridlines on the *x* and *y* axis.
        xgridminor, ygridminor : bool, optional
            Whether to draw minor gridlines for the *x* and *y* axis.
        xticklabeldir, yticklabeldir : {'in', 'out'}
            Whether to place *x* and *y* axis tick label text inside
            or outside the axes.
        xlocator, ylocator : locator spec, optional
            Used to determine the *x* and *y* axis tick mark positions. Passed
            to the `~proplot.axistools.Locator` constructor.
        xticks, yticks : optional
            Aliases for `xlocator`, `ylocator`.
        xlocator_kw, ylocator_kw : dict-like, optional
            The *x* and *y* axis locator settings. Passed to
            `~proplot.axistools.Locator`.
        xminorlocator, yminorlocator : optional
            As for `xlocator`, `ylocator`, but for the minor ticks.
        xminorticks, yminorticks : optional
            Aliases for `xminorlocator`, `yminorlocator`.
        xminorlocator_kw, yminorlocator_kw
            As for `xlocator_kw`, `ylocator_kw`, but for the minor locator.
        xformatter, yformatter : formatter spec, optional
            Used to determine the *x* and *y* axis tick label string format.
            Passed to the `~proplot.axistools.Formatter` constructor.
            Use ``[]`` or ``'null'`` for no ticks.
        xticklabels, yticklabels : optional
            Aliases for `xformatter`, `yformatter`.
        xformatter_kw, yformatter_kw : dict-like, optional
            The *x* and *y* axis formatter settings. Passed to
            `~proplot.axistools.Formatter`.
        xrotation, yrotation : float, optional
            The rotation for *x* and *y* axis tick labels. Defaults to ``0``
            for normal axes, ``rc['axes.formatter.timerotation']`` for time *x* axes.
        xtickrange, ytickrange : (float, float), optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``(-1,1)`` with
            axis range ``(-5,5)`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1.
        xmargin, ymargin : float, optional
            The default margin between plotted content and the *x* and *y*
            axis spines. Value is proportional to the width, height of the axes.
            Use this if you want whitespace between plotted content
            and the spines, but don't want to explicitly set `xlim` or `ylim`.
        xbounds, ybounds : (float, float), optional
            The *x* and *y* axis data bounds within which to draw the spines.
            For example, the axis range ``(0, 4)`` with bounds ``(1, 4)``
            will prevent the spines from meeting at the origin.
        xcolor, ycolor : color-spec, optional
            Color for the *x* and *y* axis spines, ticks, tick labels, and axis
            labels. Defaults to ``rc['color']``. Use e.g. ``ax.format(color='red')``
            to set for both axes.
        xticklen, yticklen : float or str, optional
            Tick lengths for the *x* and *y* axis. If float, units are points.
            If string, units are interpreted by `~proplot.utils.units`. Defaults
            to ``rc['ticklen']``. Minor tick lengths are scaled according
            to ``rc['ticklenratio']``. Use e.g. ``ax.format(ticklen=1)`` to
            set for both axes.
        fixticks : bool, optional
            Whether to always transform the tick locators to a
            `~matplotlib.ticker.FixedLocator` instance. Defaults to ``False``.
            If your axis ticks are doing weird things (for example, ticks
            drawn outside of the axis spine), try setting this to ``True``.
        patch_kw : dict-like, optional
            Keyword arguments used to update the background patch object. You
            can use this, for example, to set background hatching with
            ``patch_kw={'hatch':'xxx'}``.
        **kwargs
            Passed to `Axes.format` and `Axes.context`.

        Note
        ----
        If you plot something with a `numpy`
        `datetime64 <https://docs.scipy.org/doc/numpy/reference/arrays.datetime.html>`__,
        `pandas.Timestamp`, `pandas.DatetimeIndex`, `datetime.date`,
        `datetime.time`, or `datetime.datetime` array as the *x* or *y*-axis
        coordinate, the axis ticks and tick labels will be automatically
        formatted as dates.

        See also
        --------
        `~proplot.axistools.Scale`, `~proplot.axistools.Locator`,
        `~proplot.axistools.Formatter`
        """
        context, kwargs = self.context(**kwargs)
        with context:
            # Background basics
            self.patch.set_clip_on(False)
            self.patch.set_zorder(-1)
            kw_face = rc.fill({
                'facecolor': 'axes.facecolor',
                'alpha': 'axes.alpha'
                })
            patch_kw = patch_kw or {}
            kw_face.update(patch_kw)
            self.patch.update(kw_face)

            # No mutable default args
            xlabel_kw        = xlabel_kw or {}
            ylabel_kw        = ylabel_kw or {}
            xscale_kw        = xscale_kw or {}
            yscale_kw        = yscale_kw or {}
            xlocator_kw      = xlocator_kw or {}
            ylocator_kw      = ylocator_kw or {}
            xformatter_kw    = xformatter_kw or {}
            yformatter_kw    = yformatter_kw or {}
            xminorlocator_kw = xminorlocator_kw or {}
            yminorlocator_kw = yminorlocator_kw or {}
            # Flexible keyword args, declare defaults
            xmargin       = _notNone(xmargin, rc['axes.xmargin'])
            ymargin       = _notNone(ymargin, rc['axes.ymargin'])
            xtickdir      = _notNone(xtickdir, rc['xtick.direction'])
            ytickdir      = _notNone(ytickdir, rc['ytick.direction'])
            xtickminor    = _notNone(xtickminor, rc['xtick.minor.visible'])
            ytickminor    = _notNone(ytickminor, rc['ytick.minor.visible'])
            xformatter    = _notNone(xticklabels, xformatter, None, names=('xticklabels', 'xformatter'))
            yformatter    = _notNone(yticklabels, yformatter, None, names=('yticklabels', 'yformatter'))
            xlocator      = _notNone(xticks, xlocator, None, names=('xticks', 'xlocator'))
            ylocator      = _notNone(yticks, ylocator, None, names=('yticks', 'ylocator'))
            xminorlocator = _notNone(xminorticks, xminorlocator, None, names=('xminorticks', 'xminorlocator'))
            yminorlocator = _notNone(yminorticks, yminorlocator, None, names=('yminorticks', 'yminorlocator'))
            # Grid defaults are more complicated
            axis = rc.get('axes.grid.axis') # always need this property
            grid, which = rc['axes.grid'], rc['axes.grid.which']
            if which is not None or grid is not None: # if *one* was changed
                if grid is None:
                    grid = rc.get('axes.grid')
                elif which is None:
                    which = rc.get('axes.grid.which')
                xgrid = _notNone(xgrid, grid
                    and axis in ('x','both') and which in ('major','both'))
                ygrid = _notNone(ygrid, grid
                    and axis in ('y','both') and which in ('major','both'))
                xgridminor = _notNone(xgridminor, grid
                    and axis in ('x','both') and which in ('minor','both'))
                ygridminor = _notNone(ygridminor, grid
                    and axis in ('y','both') and which in ('minor','both'))

            # Sensible defaults for spine, tick, tick label, and label locs
            # NOTE: Allow tick labels to be present without ticks! User may
            # want this sometimes! Same goes for spines!
            xspineloc  = _notNone(xloc, xspineloc, None, names=('xloc', 'xspineloc'))
            yspineloc  = _notNone(yloc, yspineloc, None, names=('yloc', 'yspineloc'))
            xtickloc   = _notNone(xtickloc, xspineloc, _rcloc_to_stringloc('x', 'xtick'))
            ytickloc   = _notNone(ytickloc, yspineloc, _rcloc_to_stringloc('y', 'ytick'))
            xspineloc  = _notNone(xspineloc, _rcloc_to_stringloc('x', 'axes.spines'))
            yspineloc  = _notNone(yspineloc, _rcloc_to_stringloc('y', 'axes.spines'))
            if xtickloc != 'both':
                xticklabelloc = _notNone(xticklabelloc, xtickloc)
                xlabelloc     = _notNone(xlabelloc, xticklabelloc)
                if xlabelloc not in (None,'bottom','top'): # "both", "neither", etc
                    xlabelloc = 'bottom'
            if ytickloc != 'both':
                yticklabelloc = _notNone(yticklabelloc, ytickloc)
                ylabelloc     = _notNone(ylabelloc, yticklabelloc)
                if ylabelloc not in (None,'left','right'):
                    ylabelloc = 'left'

            # Begin loop
            for (x, axis,
                label, color, ticklen,
                margin, bounds,
                tickloc, spineloc,
                ticklabelloc, labelloc,
                grid, gridminor,
                tickminor, tickminorlocator,
                lim, reverse, scale,
                locator, tickrange,
                formatter, tickdir,
                ticklabeldir, rotation,
                label_kw, scale_kw,
                locator_kw, minorlocator_kw,
                formatter_kw
                ) in zip(
                ('x','y'), (self.xaxis, self.yaxis),
                (xlabel, ylabel), (xcolor, ycolor), (xticklen, yticklen),
                (xmargin, ymargin), (xbounds, ybounds),
                (xtickloc, ytickloc), (xspineloc, yspineloc),
                (xticklabelloc, yticklabelloc), (xlabelloc, ylabelloc),
                (xgrid, ygrid), (xgridminor, ygridminor),
                (xtickminor, ytickminor), (xminorlocator, yminorlocator),
                (xlim, ylim), (xreverse, yreverse), (xscale, yscale),
                (xlocator, ylocator), (xtickrange, ytickrange),
                (xformatter, yformatter), (xtickdir, ytickdir),
                (xticklabeldir, yticklabeldir), (xrotation, yrotation),
                (xlabel_kw, ylabel_kw), (xscale_kw, yscale_kw),
                (xlocator_kw, ylocator_kw), (xminorlocator_kw, yminorlocator_kw),
                (xformatter_kw, yformatter_kw),
                ):
                # Axis scale
                # WARNING: This relies on monkey patch of mscale.scale_factory
                # that allows it to accept a custom scale class!
                if scale is not None:
                    if (formatter is None and getattr(scale,'name',scale) in
                        ('log','logit','inverse','symlog')):
                        formatter = 'simple'
                    scale, args, kw = axistools.Scale(scale, **scale_kw)
                    getattr(self, f'set_{x}scale')(scale, *args, **kw)
                # Axis limits
                # NOTE: 3.1+ has axis.set_inverted(), below is from source code
                if lim is not None:
                    getattr(self, f'set_{x}lim')(lim)
                if reverse:
                    lo, hi = axis.get_view_interval()
                    axis.set_view_interval(
                        max(lo, hi), min(lo, hi), ignore=True)
                # Is this a date axis?
                # NOTE: Make sure to get this *after* lims set!
                date = isinstance(axis.converter, mdates.DateConverter)

                # Fix spines
                kw = rc.fill({
                    'linewidth': 'axes.linewidth',
                    'color':     'axes.edgecolor',
                    })
                if color is not None:
                    kw['color'] = color
                sides = ('bottom','top') if x == 'x' else ('left','right')
                spines = [self.spines[s] for s in sides]
                for spine,side in zip(spines,sides):
                    # Line properties
                    # Override if we're settings spine bounds
                    if bounds is not None and spineloc not in sides:
                        spineloc = sides[0] # by default, should just have spines on edges in this case
                    # Eliminate sides
                    if spineloc == 'neither':
                        spine.set_visible(False)
                    elif spineloc == 'both':
                        spine.set_visible(True)
                    elif spineloc in sides: # make relevant spine visible
                        b = True if side == spineloc else False
                        spine.set_visible(b)
                    elif spineloc is not None:
                        # Special spine location, usually 'zero', 'center',
                        # or tuple with (units, location) where 'units' can
                        # be 'axes', 'data', or 'outward'.
                        if side == sides[1]:
                            spine.set_visible(False)
                        else:
                            spine.set_visible(True)
                            try:
                                spine.set_position(spineloc)
                            except ValueError:
                                raise ValueError(f'Invalid {x} spine location {spineloc!r}. Options are {", ".join((*sides, "both", "neither"))}.')
                    # Apply spine bounds
                    if bounds is not None and spine.get_visible():
                        spine.set_bounds(*bounds)
                    spine.update(kw)
                # Get which spines are visible; needed for setting tick locations
                spines = [side for side,spine in zip(sides,spines) if spine.get_visible()]

                # Tick and grid settings for major and minor ticks separately
                # Override is just a "new default", but user can override this
                grid_dict = lambda grid: {
                    'grid_color':     grid + '.color',
                    'grid_alpha':     grid + '.alpha',
                    'grid_linewidth': grid + '.linewidth',
                    'grid_linestyle': grid + '.linestyle',
                    }
                for which,igrid in zip(('major', 'minor'), (grid,gridminor)):
                    # Tick properties
                    kw_ticks = rc.category(x + 'tick.' + which)
                    if kw_ticks is None:
                        kw_ticks = {}
                    else:
                        kw_ticks.pop('visible', None) # invalid setting
                    if ticklen is not None:
                        if which == 'major':
                            kw_ticks['size'] = utils.units(ticklen, 'pt')
                        else:
                            kw_ticks['size'] = utils.units(ticklen, 'pt') * rc.get('ticklenratio')
                    # Grid style and toggling
                    if igrid is not None:
                        axis.grid(igrid, which=which) # toggle with special global props
                    if which == 'major':
                        kw_grid = rc.fill(grid_dict('grid'))
                    else:
                        kw_major = kw_grid
                        kw_grid = rc.fill(grid_dict('gridminor'))
                        kw_grid.update({key:value for key,value in kw_major.items() if key not in kw_grid})
                    # Changed rc settings
                    axis.set_tick_params(which=which, **kw_grid, **kw_ticks)

                # Tick and ticklabel properties that apply to both major and minor
                # * Weird issue seems to cause set_tick_params to reset/forget that the grid
                #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
                #   Since gridOn is undocumented feature, don't use it. So calling _format_axes() a second time will remove the lines
                # * Can specify whether the left/right/bottom/top spines get ticks; sides will be
                #   group of left/right or top/bottom
                # * Includes option to draw spines but not draw ticks on that spine, e.g.
                #   on the left/right edges
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
                        if len(options) == 1:
                            labelloc = options[0]
                elif labelloc not in sides:
                    raise ValueError(f'Got labelloc {labelloc!r}, valid options are {sides!r}.')
                # Apply
                axis.set_tick_params(which='both', **kw)
                if labelloc is not None:
                    axis.set_label_position(labelloc)

                # Tick label settings
                # First color and size
                kw = rc.fill({
                    'labelcolor': 'tick.labelcolor', # new props
                    'labelsize': 'tick.labelsize',
                    'color': x + 'tick.color',
                    })
                if color:
                    kw['color'] = color
                    kw['labelcolor'] = color
                # Tick direction and rotation
                if tickdir == 'in':
                    kw['pad'] = 1 # ticklabels should be much closer
                if ticklabeldir == 'in': # put tick labels inside the plot
                    tickdir = 'in'
                    pad = (rc.get(x + 'tick.major.size')
                           + rc.get(x + 'tick.major.pad')
                           + rc.get(x + 'tick.labelsize'))
                    kw['pad'] = -pad
                if tickdir is not None:
                    kw['direction'] = tickdir
                axis.set_tick_params(which='both', **kw)

                # Settings that can't be controlled by set_tick_params
                # Also set rotation and alignment here
                kw = rc.fill({
                    'fontfamily': 'font.family',
                    'weight':     'tick.labelweight'
                    })
                if rotation is not None:
                    kw = {'rotation':rotation}
                    if x == 'x':
                        self._datex_rotated = True
                        if rotation not in (0,90,-90):
                            kw['ha'] = ('right' if rotation > 0 else 'left')
                for t in axis.get_ticklabels():
                    t.update(kw)
                # Margins
                if margin is not None:
                    self.margins(**{x: margin})

                # Axis label updates
                # NOTE: This has to come after set_label_position, or ha or va
                # overrides in label_kw are overwritten
                kw = rc.fill({
                    'color':      'axes.edgecolor',
                    'fontsize':   'axes.labelsize',
                    'weight':     'axes.labelweight',
                    'fontfamily': 'font.family',
                    })
                if label is not None:
                    kw['text'] = label
                if color:
                    kw['color'] = color
                kw.update(label_kw)
                if kw: # NOTE: initially keep spanning labels off
                    self.figure._update_axislabels(axis, **kw)

                # Major and minor locator
                # WARNING: MultipleLocator fails sometimes, notably when doing
                # boxplot. Tick labels moved to left and are incorrect.
                if locator is not None:
                    locator = axistools.Locator(locator, **locator_kw)
                    axis.set_major_locator(locator)
                    if isinstance(locator, mticker.IndexLocator):
                        tickminor = False # minor ticks make no sense for 'index' data
                if not tickminor and tickminorlocator is None:
                    axis.set_minor_locator(axistools.Locator('null'))
                elif tickminorlocator is not None:
                    axis.set_minor_locator(axistools.Locator(tickminorlocator, **minorlocator_kw))

                # Major and minor formatter
                # NOTE: Only reliable way to disable ticks labels and then
                # restore them is by messing with the formatter, *not* setting
                # labelleft=False, labelright=False, etc. Check for this here
                fixedformatfix = False
                if formatter is not None or tickrange is not None and not (
                    isinstance(axis.get_major_formatter(), mticker.NullFormatter)
                    and getattr(self, '_share' + x)):
                    # Tick range
                    if tickrange is not None:
                        if formatter not in (None,'auto'):
                            warnings.warn('The tickrange feature requires proplot.AutoFormatter formatter. Overriding input formatter.')
                        formatter = 'auto'
                        formatter_kw.setdefault('tickrange', tickrange)
                    # Set the formatter
                    if formatter in ('date','concise'):
                        locator = axis.get_major_locator()
                        formatter_kw.setdefault('locator', locator)
                    if (isinstance(axis.get_major_formatter(), mticker.NullFormatter)
                        and getattr(self, '_share' + x)):
                        pass # this is a shared axis with disabled ticks
                    else:
                        formatter = axistools.Formatter(formatter, date=date, **formatter_kw)
                        axis.set_major_formatter(formatter)
                    if isinstance(formatter, mticker.FixedFormatter): # if locator is MultipleLocator, first tick gets cut off!
                        fixedformatfix = True
                axis.set_minor_formatter(mticker.NullFormatter())

                # Ensure no out-of-bounds ticks! Even set_smart_bounds()
                # fails sometimes.
                # * Using set_bounds also failed, and fancy method overrides did
                #   not work, so instead just turn locators into fixed version
                # * Most locators take no arguments in __call__, and some do not
                #   have tick_values method, so we just call them.
                if fixticks or fixedformatfix or bounds is not None or axis.get_scale() == 'cutoff':
                    if bounds is None:
                        bounds = getattr(self, f'get_{x}lim')()
                    locator = axistools.Locator([x for x in axis.get_major_locator()() if bounds[0] <= x <= bounds[1]])
                    axis.set_major_locator(locator)
                    locator = axistools.Locator([x for x in axis.get_minor_locator()() if bounds[0] <= x <= bounds[1]])
                    axis.set_minor_locator(locator)
            # Call parent
            if aspect is not None:
                self.set_aspect(aspect)
            super().format(**kwargs)

    def altx(self, *args, **kwargs):
        """Alias (and more intuitive name) for `~CartesianAxes.twiny`.
        The matplotlib `~matplotlib.axes.Axes.twiny` function
        generates two *x*-axes with a shared ("twin") *y*-axis."""
        # Cannot wrap twiny() because we want to use CartesianAxes, not
        # matplotlib Axes. Instead use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py
        if self._altx_child:
            raise ValueError('No more than *two* twin axes!')
        if self._altx_parent:
            raise ValueError('This *is* a twin axes!')
        with self.figure._unlock():
            ax = self._make_twin_axes(sharey=self, projection='cartesian')
        # ax.set_autoscaley_on(self.get_autoscaley_on()) # shared axes must have matching autoscale
        ax.grid(False)
        self._altx_child = ax
        ax._altx_parent = self
        self._altx_overrides()
        ax._altx_overrides()
        self.add_child_axes(ax)
        return ax

    def alty(self):
        """Alias (and more intuitive name) for `~CartesianAxes.twinx`.
        The matplotlib `~matplotlib.axes.Axes.twinx` function
        generates two *y*-axes with a shared ("twin") *x*-axis."""
        # Must reproduce twinx here because need to generate CartesianAxes
        if self._alty_child:
            raise ValueError('No more than *two* twin axes!')
        if self._alty_parent:
            raise ValueError('This *is* a twin axes!')
        with self.figure._unlock():
            ax = self._make_twin_axes(sharex=self, projection='cartesian')
        # ax.set_autoscalex_on(self.get_autoscalex_on()) # shared axes must have matching autoscale
        ax.grid(False)
        self._alty_child = ax
        ax._alty_parent = self
        self._alty_overrides()
        ax._alty_overrides()
        self.add_child_axes(ax)
        return ax

    def dualx(self, offset=0, scale=1, xscale='linear', xlabel=None, **kwargs):
        """
        Makes a secondary *x*-axis for denoting equivalent *x*
        coordinates in *alternate units*. Returns nothing.

        Parameters
        ----------
        scale : float, optional
            The constant multiple applied after scaling data with `transform`.
            Defaults to ``1``.
            For example, if your *x*-axis is meters and you
            want kilometers on the other side, use ``scale=1e-3``.
        offset : float, optional
            The constant offset added after multipyling by `scale`.
            Defaults to ``0``.
            For example, if your *x*-axis is Kelvin and you want degrees
            Celsius on the opposite side, use ``offset=-273.15``.
        xscale : str, optional
            The registered scale name used to transform data to the alternate
            units.  Defaults to ``'linear'``.
            For example, if your *x*-axis is wavenumber and you want wavelength
            on the opposite side, use ``xscale='inverse'``. If your *x*-axis
            is height and you want pressure on the opposite side, use
            ``xscale='pressure'``. For the opposite, use ``xscale='height'``.
        xlabel : str, optional
            The axis label (highly recommended). A warning will be issued if
            this is not supplied.
        **kwargs
            Passed to `Axes.format`.
        """
        # The axis scale is used to transform units on the left axis, linearly
        # spaced, to units on the right axis... so the right scale must scale
        # its data with the *inverse* of this transform. We do this below.
        # TODO: Use matplotlib 3.1 'secondary axis' feature? Adds a 'child axes'
        # just like 'insets' -- also *colorbars* (these are only two ways that
        # child axes are added). *Only* place child axes appear is in automatic
        # title positioning and in get_children function.
        parent = self.get_xscale()
        if parent != 'linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current {parent!r} scale.')
            self.set_xscale('linear')
        ax = self.twiny()
        if xlabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "xlabel" keyword argument.')
        xscale = axistools.InvertedScaleFactory(xscale)
        ax.format(xscale=xscale, xlabel=xlabel, **kwargs)
        self._dualx_scale = (offset, scale)
        return ax

    def dualy(self, offset=0, scale=1, yscale='linear', ylabel=None, **kwargs):
        """
        Makes a secondary *y*-axis for denoting equivalent *y*
        coordinates in *alternate units*. Returns nothing.

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
            For example, if your *y*-axis is wavenumber and you want wavelength
            on the opposite side, use ``yscale='inverse'``. If your *y*-axis
            is height and you want pressure on the opposite side, use
            ``yscale='pressure'``. For the opposite, use ``xscale='height'``.
        ylabel : str, optional
            The axis label (highly recommended). A warning will be issued if
            this is not supplied.
        **kwargs
            Passed to `Axes.format`.
        """
        parent = self.get_yscale()
        if parent != 'linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current {parent!r} scale.')
            self.set_yscale('linear')
        ax = self.twinx()
        if ylabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "ylabel" keyword argument.')
        yscale = axistools.InvertedScaleFactory(yscale)
        ax.format(yscale=yscale, ylabel=ylabel, **kwargs)
        self._dualy_scale = (offset, scale)
        return ax

    def draw(self, renderer=None, *args, **kwargs):
        """Adds post-processing steps before axes is drawn."""
        # NOTE: This mimics matplotlib API, which calls identical
        # post-processing steps in both draw() and get_tightbbox()
        self._datex_rotate()
        self._dualx_lock()
        self._dualy_lock()
        self._altx_overrides()
        self._alty_overrides()
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        super().draw(renderer, *args, **kwargs)

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Adds post-processing steps before tight bounding box is
        calculated."""
        self._datex_rotate()
        self._dualx_lock()
        self._dualy_lock()
        self._altx_overrides()
        self._alty_overrides()
        return super().get_tightbbox(renderer, *args, **kwargs)

    def twinx(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twinx` and intelligently
        handles axis ticks, gridlines, axis tick labels, axis labels, and axis
        sharing. Returns a `CartesianAxes` instance."""
        return self.alty()

    def twiny(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twiny` and intelligently
        handles axis ticks, gridlines, axis tick labels, axis labels, and axis
        sharing. Returns a `CartesianAxes` instance."""
        return self.altx()

class ProjectionAxes(Axes):
    """Intermediate class, shared by `CartopyAxes` and
    `BasemapAxes`. Disables methods that are inappropriate for map
    projections and adds `ProjectionAxes.format`, so that arguments
    passed to `Axes.format` are identical for `CartopyAxes`
    and `BasemapAxes`."""
    def __init__(self, *args, **kwargs): # just to disable docstring inheritence
        """
        See also
        --------
        `~proplot.subplots.subplots`, `CartopyAxes`, `BasemapAxes`
        """
        super().__init__(*args, **kwargs)

    @wrappers._expand_methods_list
    def __getattribute__(self, attr, *args):
        """Disables the methods `MAP_DISABLED_METHODS`, which are
        inappropriate for map projections."""
        if attr in wrappers.MAP_DISABLED_METHODS:
            raise AttributeError(f'Invalid plotting function {attr!r} for map projection axes.')
        return super().__getattribute__(attr, *args)

    def _projection_format_kwargs(self, *,
        labels=None, latlabels=None, lonlabels=None,
        lonlim=None, latlim=None, latmax=None, grid=None,
        lonlocator=None, lonlines=None,
        latlocator=None, latlines=None,
        boundinglat=None,
        **kwargs,
        ):
        # Parse alternative keyword args
        latmax = _notNone(latmax, rc['geogrid.latmax'])
        lonlines = _notNone(lonlines, lonlocator, rc['geogrid.lonstep'], names=('lonlines', 'lonlocator'))
        latlines = _notNone(latlines, latlocator, rc['geogrid.latstep'], names=('latlines', 'latlocator'))
        if lonlabels is not None or latlabels is not None:
            labels = True
        else:
            labels = _notNone(labels, rc['geogrid.labels'])
        if (lonlines is not None or latlines is not None or
                latmax is not None or labels):
            grid = True
        else:
            grid = _notNone(grid, rc['geogrid'])

        # Longitude gridlines, draw relative to projection prime meridian
        if isinstance(self, CartopyAxes):
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
        else:
            base = 5
            lon_0 = base*round(self.projection.lonmin/base) + 180 # central longitude
        if lonlines is not None:
            if not np.iterable(lonlines):
                lonlines = utils.arange(lon_0 - 180, lon_0 + 180, lonlines)
            lonlines = [*lonlines]
        # Latitudes gridlines, draw from -latmax to latmax, but if result would
        # be asymmetrical across equator, do not use
        if latlines is not None or latmax is not None:
            # Fill defaults
            if latlines is None:
                latlines = _notNone(self._latlines_values, rc.get('geogrid.latstep'))
            if latmax is None:
                latmax = _notNone(self._latmax, rc.get('geogrid.latmax'))
            # Get tick locations
            if not np.iterable(latlines):
                if (latmax % latlines) == (-latmax % latlines):
                    latlines = utils.arange(-latmax, latmax, latlines)
                else:
                    latlines = utils.arange(0, latmax, latlines)
                    if latlines[-1] != latmax:
                        latlines = np.concatenate((latlines, [latmax]))
                    latlines = np.concatenate((-latlines[::-1], latlines[1:]))
            latlines = [*latlines]
        # Add attributes
        if latmax is not None:
            self._latmax = latmax
        if latlines is not None:
            self._latlines_values = latlines
        if lonlines is not None:
            self._lonlines_values = lonlines

        # Length-4 boolean arrays of whether and where to toggle labels
        # Format is [left, right, bottom, top]
        if lonlabels or latlabels:
            labels = True # toggle them at all?
        lonarray, latarray = [], [] # destination
        for labs,array in zip((lonlabels,latlabels), (lonarray,latarray)):
            if labs is False:
                return [0]*4
            if labs is None:
                labs = 1
            if isinstance(labs, str):
                string = labs
                labs = [0]*4
                for idx,char in zip([0,1,2,3],'lrbt'):
                    if char in string:
                        labs[idx] = 1
            if isinstance(labs, Number): # e.g. *boolean*
                labs = np.atleast_1d(labs)
            if len(labs) == 1:
                labs = [*labs, 0] # default is to label bottom/left
            if len(labs) == 2:
                if array is lonarray:
                    labs = [0, 0, *labs]
                else:
                    labs = [*labs, 0, 0]
            elif len(labs) != 4:
                raise ValueError(f'Invalid lon/lat label spec: {labs}.')
            array[:] = labs
        return (grid, latmax, lonlim, latlim, boundinglat,
                lonlines, latlines, labels, lonarray, latarray, kwargs)

class PolarAxes(ProjectionAxes, mproj.PolarAxes):
    """Intermediate class, mixes `ProjectionAxes` with
    `~matplotlib.projections.polar.PolarAxes`."""
    name = 'polar2'
    """The registered projection name."""
    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`
        """
        # Set tick length to zero so azimuthal labels are not too offset
        # Change default radial axis formatter but keep default theta one
        super().__init__(*args, **kwargs)
        formatter = axistools.Formatter('default')
        self.yaxis.set_major_formatter(formatter)
        self.yaxis.isDefault_majfmt = True
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)

    def __getattribute__(self, attr, *args):
        """Applies the `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`,
        `~proplot.wrappers.enforce_centers`, `~proplot.wrappers.enforce_edges`,
        `~proplot.wrappers.cartopy_gridfix`, `~proplot.wrappers.cartopy_transform`,
        `~proplot.wrappers.cartopy_crs`, `~proplot.wrappers.plot_wrapper`,
        `~proplot.wrappers.scatter_wrapper`, `~proplot.wrappers.fill_between_wrapper`,
        and `~proplot.wrappers.fill_betweenx_wrapper` wrappers."""
        obj = super().__getattribute__(attr, *args)
        if callable(obj):
            # Step 4) Color usage wrappers
            if attr in wrappers.CMAP_METHODS:
                obj = wrappers._cmap_wrapper(self, obj)
            elif attr in wrappers.CYCLE_METHODS:
                obj = wrappers._cycle_wrapper(self, obj)
            # Step 3) Fix coordinate grid
            if attr in wrappers.EDGES_METHODS or attr in wrappers.CENTERS_METHODS:
                obj = wrappers._cartopy_gridfix(self, obj)
            # Step 2) Utilities
            if attr in wrappers.EDGES_METHODS:
                obj = wrappers._enforce_edges(self, obj)
            elif attr in wrappers.CENTERS_METHODS:
                obj = wrappers._enforce_centers(self, obj)
            # Step 1) Parse args input
            if attr in wrappers.D2_METHODS:
                obj = wrappers._autoformat_2d(self, obj)
            elif attr in wrappers.D1_METHODS:
                obj = wrappers._autoformat_1d(self, obj)
            # Step 0) Special wrappers
            if attr == 'plot':
                obj = wrappers._plot_wrapper(self, obj)
            elif attr == 'scatter':
                obj = wrappers._scatter_wrapper(self, obj)
            elif attr == 'fill_between':
                obj = wrappers._fill_between_wrapper(self, obj)
            elif attr == 'fill_betweenx':
                obj = wrappers._fill_betweenx_wrapper(self, obj)
        return obj

    def format(self, *args,
        r0=None, theta0=None, thetadir=None,
        thetamin=None, thetamax=None, thetalim=None,
        rmin=None, rmax=None, rlim=None,
        rlabelpos=None, rscale=None, rborder=None,
        thetalocator=None, rlocator=None, thetalines=None, rlines=None,
        thetaformatter=None, rformatter=None, thetalabels=None, rlabels=None,
        thetalocator_kw=None, rlocator_kw=None,
        thetaformatter_kw=None, rformatter_kw=None,
        **kwargs):
        """
        Calls `CartesianAxes.format` and formats the tick locations, tick
        labels, grid lines, and more. All ``theta`` arguments are
        specified in **degrees**, not radians. The below parameters are
        specific to `PolarAxes`.

        Parameters
        ----------
        r0 : float, optional
            The radial origin.
        theta0 : {'N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'}
            The zero azimuth location.
        thetadir : {-1, 1, 'clockwise', 'anticlockwise', 'counterclockwise'}, optional
            The positive azimuth direction. Clockwise corresponds to ``-1``
            and anticlockwise corresponds to ``-1``. Defaults to ``-1``.
        thetamin, thetamax : float, optional
            The lower and upper azimuthal bounds in degrees. If
            ``thetamax != thetamin + 360``, this produces a sector plot.
        thetalim : (float, float), optional
            Specifies `thetamin` and `thetamax` at once.
        rmin, rmax : float, optional
            The inner and outer radial limits. If ``r0 != rmin``, this
            produces an annular plot.
        rlim : (float, float), optional
            Specifies `rmin` and `rmax` at once.
        rborder : bool, optional
            Toggles the polar axes border on and off. Visibility of the "inner"
            radial spine and "start" and "end" azimuthal spines is controlled
            automatically be matplotlib.
        thetalocator, rlocator : float or list of float, optional
            Used to determine the azimuthal and radial gridline positions.
            Passed to the `~proplot.axistools.Locator` constructor.
        thetalines, rlines
            Aliases for `thetalocator`, `rlocator`.
        thetalocator_kw, rlocator_kw : dict-like, optional
            The azimuthal and radial locator settings. Passed to
            `~proplot.axistools.Locator`.
        rlabelpos : float, optional
            The azimuth at which radial coordinates are labeled.
        thetaformatter, rformatter : formatter spec, optional
            Used to determine the azimuthal and radial label format.
            Passed to the `~proplot.axistools.Formatter` constructor.
            Use ``[]`` or ``'null'`` for no ticks.
        thetalabels, rlabels : optional
            Aliases for `thetaformatter`, `rformatter`.
        thetaformatter_kw, rformatter_kw : dict-like, optional
            The azimuthal and radial label formatter settings. Passed to
            `~proplot.axistools.Formatter`.
        **kwargs
            Passed to `Axes.format` and `Axes.context`
        """
        context, kwargs = self.context(**kwargs)
        with context:
            # Not mutable default args
            thetalocator_kw   = thetalocator_kw or {}
            thetaformatter_kw = thetaformatter_kw or {}
            rlocator_kw       = rlocator_kw or {}
            rformatter_kw     = rformatter_kw or {}
            # Flexible input
            if rlim is not None:
                if rmin is not None or rmax is not None:
                    warnings.warn(f'Conflicting keyword args rmin={rmin}, rmax={rmax}, and rlim={rlim}. Using "rlim".')
                rmin, rmax = rlim
            if thetalim is not None:
                if thetamin is not None or thetamax is not None:
                    warnings.warn(f'Conflicting keyword args thetamin={thetamin}, thetamax={thetamax}, and thetalim={thetalim}. Using "thetalim".')
                thetamin, thetamax = thetalim
            thetalocator   = _notNone(thetalines, thetalocator, None, names=('thetalines', 'thetalocator'))
            thetaformatter = _notNone(thetalabels, thetaformatter, None, names=('thetalabels', 'thetaformatter'))
            rlocator       = _notNone(rlines, rlocator, None, names=('rlines', 'rlocator'))
            rformatter     = _notNone(rlabels, rformatter, None, names=('rlabels', 'rformatter'))

            # Special radius settings
            if r0 is not None:
                self.set_rorigin(r0)
            if rlabelpos is not None:
                self.set_rlabel_position(rlabelpos)
            if rscale is not None:
                self.set_rscale(rscale)
            if rborder is not None:
                self.spines['polar'].set_visible(bool(rborder))
            # Special azimuth settings
            if theta0 is not None:
                self.set_theta_zero_location(theta0)
            if thetadir is not None:
                self.set_theta_direction(thetadir)

            # Iterate
            for (x, name, axis,
                min_, max_,
                locator, formatter,
                locator_kw, formatter_kw,
                ) in zip(
                ('x','y'), ('theta', 'r'), (self.xaxis, self.yaxis),
                (thetamin, rmin), (thetamax, rmax),
                (thetalocator, rlocator), (thetaformatter, rformatter),
                (thetalocator_kw, rlocator_kw), (thetaformatter_kw, rformatter_kw),
                ):
                # Axis limits
                # Try to use public API where possible
                if min_ is not None:
                    getattr(self, 'set_' + name + 'min')(min_)
                else:
                    min_ = getattr(self, 'get_' + name + 'min')()
                if max_ is not None:
                    getattr(self, 'set_' + name + 'max')(max_)
                else:
                    max_ = getattr(self, 'get_' + name + 'max')()

                # Spine settings
                kw = rc.fill({
                    'linewidth': 'axes.linewidth',
                    'color': 'axes.edgecolor',
                    })
                sides = ('inner','polar') if name == 'r' else ('start','end')
                spines = [self.spines[s] for s in sides]
                for spine,side in zip(spines,sides):
                    spine.update(kw)

                # Grid and grid label settings
                # NOTE: Not sure if polar lines inherit tick or grid props
                kw = rc.fill({
                    'color': x + 'tick.color',
                    'labelcolor': 'tick.labelcolor', # new props
                    'labelsize': 'tick.labelsize',
                    'grid_color': 'grid.color',
                    'grid_alpha': 'grid.alpha',
                    'grid_linewidth': 'grid.linewidth',
                    'grid_linestyle': 'grid.linestyle',
                    })
                axis.set_tick_params(which='both', **kw)
                # Label settings that can't be controlled with set_tick_params
                kw = rc.fill({
                    'fontfamily': 'font.family',
                    'weight':     'tick.labelweight'
                    })
                for t in axis.get_ticklabels():
                    t.update(kw)

                # Tick locator, which in this case applies to gridlines
                # NOTE: Must convert theta locator input to radians, then back
                # to degrees.
                if locator is not None:
                    if name == 'theta' and (
                        not isinstance(locator, (str,mticker.Locator))):
                        locator = np.deg2rad(locator) # real axis limts are rad
                    locator = axistools.Locator(locator, **locator_kw)
                    locator.set_axis(axis) # this is what set_locator does
                    grids = np.array(locator())
                    if name == 'r':
                        grids = grids[(grids >= min_) & (grids <= max_)]
                        self.set_rgrids(grids)
                    else:
                        grids = np.rad2deg(grids)
                        grids = grids[(grids >= min_) & (grids <= max_)]
                        if grids[-1] == min_+360: # exclusive if 360 degrees
                            grids = grids[:-1]
                        self.set_thetagrids(grids)
                # Tick formatter and toggling
                if formatter is not None:
                    formatter = axistools.Formatter(formatter, **formatter_kw)
                    axis.set_major_formatter(formatter)

            # Parent method
            super().format(*args, **kwargs)

# Cartopy takes advantage of documented feature where any class with method
# named _as_mpl_axes can be passed as 'projection' object.
# Feature documented here: https://matplotlib.org/devel/add_new_projection.html
# class CartopyAxes(ProjectionAxes, GeoAxes):
class CartopyAxes(ProjectionAxes, GeoAxes):
    """Axes subclass for plotting `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__
    projections. Initializes the `cartopy.crs.Projection` instance. Also
    allows for *partial* coverage of azimuthal projections by zooming into
    the full projection, then drawing a circle boundary around some latitude
    away from the center (this is surprisingly difficult to do)."""
    name = 'cartopy'
    """The registered projection name."""
    _n_points = 100 # number of points for drawing circle map boundary
    _proj_circles = ('laea', 'aeqd', 'stere', 'gnom')
    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        *args, **kwargs
            Passed to `Axes.__init__`.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.proj`
        """
        # GeoAxes initialization steps are run manually
        # If _hold is set False or None, cartopy will call cla() on axes,
        # which wipes out row and column labels, so we do not use it
        import cartopy.crs as ccrs
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError('You must initialize CartopyAxes with map_projection=<cartopy.crs.Projection > .')
        self.projection = map_projection # attribute used with GeoAxes
        self.img_factories = []
        self.outline_patch = None
        self.background_patch = None
        self._gridliners = [] # populated in first format call
        self._done_img_factory = False
        self._boundinglat = None # set in first format call
        self._latmax = None # used for format memory
        self._lonlines_values = None
        self._latlines_values = None
        super().__init__(*args, map_projection=map_projection, **kwargs)
        # Zero ticks so gridlines are not offset
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)
        # Default bounds and extent, and always user circle for some projs
        proj = self.projection.proj4_params['proj']
        if proj not in self._proj_circles:
            self.set_global() # see: https://stackoverflow.com/a/48956844/4970632
        else:
            self.set_boundary(projs.Circle(self._n_points),
                              transform=self.transAxes)

    def __getattribute__(self, attr, *args):
        """Applies the `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`,
        `~proplot.wrappers.enforce_centers`, `~proplot.wrappers.enforce_edges`,
        `~proplot.wrappers.cartopy_gridfix`, `~proplot.wrappers.cartopy_transform`,
        `~proplot.wrappers.cartopy_crs`, `~proplot.wrappers.plot_wrapper`,
        `~proplot.wrappers.scatter_wrapper`, `~proplot.wrappers.fill_between_wrapper`,
        and `~proplot.wrappers.fill_betweenx_wrapper` wrappers."""
        obj = super().__getattribute__(attr, *args)
        if callable(obj):
            # Step 5) Color usage wrappers
            if attr in wrappers.CMAP_METHODS:
                obj = wrappers._cmap_wrapper(self, obj)
            elif attr in wrappers.CYCLE_METHODS:
                obj = wrappers._cycle_wrapper(self, obj)
            # Step 4) Fix coordinate grid
            if attr in wrappers.EDGES_METHODS or attr in wrappers.CENTERS_METHODS:
                obj = wrappers._cartopy_gridfix(self, obj)
            # Step 3) Utilities
            if attr in wrappers.EDGES_METHODS:
                obj = wrappers._enforce_edges(self, obj)
            elif attr in wrappers.CENTERS_METHODS:
                obj = wrappers._enforce_centers(self, obj)
            # Step 2) Better default keywords
            if attr in wrappers.TRANSFORM_METHODS:
                obj = wrappers._cartopy_transform(self, obj)
            elif attr in wrappers.CRS_METHODS:
                obj = wrappers._cartopy_crs(self, obj)
            # Step 1) Parse args input
            if attr in wrappers.D2_METHODS:
                obj = wrappers._autoformat_2d(self, obj)
            elif attr in wrappers.D1_METHODS:
                obj = wrappers._autoformat_1d(self, obj)
            # Step 0) Special wrappers
            if attr == 'plot':
                obj = wrappers._plot_wrapper(self, obj)
            elif attr == 'scatter':
                obj = wrappers._scatter_wrapper(self, obj)
            elif attr == 'fill_between':
                obj = wrappers._fill_between_wrapper(self, obj)
            elif attr == 'fill_betweenx':
                obj = wrappers._fill_betweenx_wrapper(self, obj)
        return obj

    def format(self, *, patch_kw=None, **kwargs):
        # Docstring added at bottom
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs
        from cartopy.mpl import gridliner
        # Initial gridliner object, which ProPlot passively modifies
        # TODO: Flexible formatter?
        if not self._gridliners:
            gl = self.gridlines(zorder=5, draw_labels=False)
            gl.xlines = False
            gl.ylines = False
            try:
                lonformat = gridliner.LongitudeFormatter # newer
                latformat = gridliner.LatitudeFormatter
            except AttributeError:
                lonformat = gridliner.LONGITUDE_FORMATTER # older
                latformat = gridliner.LATITUDE_FORMATTER
            gl.xformatter = lonformat
            gl.yformatter = latformat
            gl.xlabels_top    = False
            gl.xlabels_bottom = False
            gl.ylabels_left   = False
            gl.ylabels_right  = False

        # Format
        context, kwargs = self.context(**kwargs)
        with context:
            (grid, _, lonlim, latlim, boundinglat,
                lonlocator, latlocator, labels, lonlabels, latlabels,
                kwargs) = self._projection_format_kwargs(**kwargs)
            # Projection extent
            # NOTE: They may add this as part of set_xlim and set_ylim in the
            # near future; see: https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638
            # WARNING: The set_extent method tries to set a *rectangle* between
            # the *4* (x,y) coordinate pairs (each corner), so something like
            # (-180,180,-90,90) will result in *line*, causing error!
            proj = self.projection.proj4_params['proj']
            north = isinstance(self.projection, (ccrs.NorthPolarStereo,
                projs.NorthPolarGnomonic,
                projs.NorthPolarAzimuthalEquidistant,
                projs.NorthPolarLambertAzimuthalEqualArea))
            south = isinstance(self.projection, (ccrs.SouthPolarStereo,
                projs.SouthPolarGnomonic,
                projs.SouthPolarAzimuthalEquidistant,
                projs.SouthPolarLambertAzimuthalEqualArea))
            if north or south:
                if (lonlim is not None or latlim is not None):
                    warnings.warn(f'{proj!r} extent is controlled by "boundinglat", ignoring lonlim={lonlim!r} and latlim={latlim!r}.')
                if self._boundinglat is None:
                    if isinstance(self.projection, projs.NorthPolarGnomonic):
                        boundinglat = 30
                    elif isinstance(self.projection, projs.SouthPolarGnomonic):
                        boundinglat = -30
                    else:
                        boundinglat = 0
                if boundinglat is not None and boundinglat != self._boundinglat:
                    eps = 1e-10 # had bug with full -180, 180 range when lon_0 not 0
                    lat0 = (90 if north else -90)
                    lon0 = self.projection.proj4_params['lon_0']
                    extent = [lon0 - 180 + eps, lon0 + 180 - eps, boundinglat, lat0]
                    self.set_extent(extent, crs=ccrs.PlateCarree())
                    self._boundinglat = boundinglat
            else:
                if boundinglat is not None:
                    warnings.warn(f'{proj!r} extent is controlled by "lonlim" and "latlim", ignoring boundinglat={boundinglat!r}.')
                if lonlim is not None or latlim is not None:
                    lonlim = lonlim or [None, None]
                    latlim = latlim or [None, None]
                    lonlim, latlim = [*lonlim], [*latlim]
                    lon_0 = self.projection.proj4_params.get('lon_0', 0)
                    if lonlim[0] is None:
                        lonlim[0] = lon_0 - 180
                    if lonlim[1] is None:
                        lonlim[1] = lon_0 + 180
                    eps = 1e-10 # had bug with full -180, 180 range when lon_0 was not 0
                    lonlim[0] += eps
                    if latlim[0] is None:
                        latlim[0] = -90
                    if latlim[1] is None:
                        latlim[1] = 90
                    extent = [*lonlim, *latlim]
                    self.set_extent(extent, crs=ccrs.PlateCarree())

            # Draw gridlines, manage them with one custom gridliner generated
            # by ProPlot, user may want to use griliner API directly
            gl = self._gridliners[0]
            if grid:
                # Collection props, see GoeAxes.gridlines() source code
                kw = rc.fill({
                    'alpha':     'geogrid.alpha',
                    'color':     'geogrid.color',
                    'linewidth': 'geogrid.linewidth',
                    'linestyle': 'geogrid.linestyle',
                    }) # cached changes
                gl.collection_kwargs.update(kw)
                # Grid locations
                # TODO: Check eps
                eps = 1e-10
                if lonlocator is not None:
                    gl.xlines = True
                    gl.xlocator = mticker.FixedLocator(lonlocator)
                if latlocator is not None:
                    gl.ylines = True
                    if latlocator[0] == -90:
                        latlocator[0] += eps
                    if latlocator[-1] == 90:
                        latlocator[-1] -= eps
                    gl.ylocator = mticker.FixedLocator(latlocator)
                if labels and not isinstance(self.projection, (ccrs.Mercator,
                                                           ccrs.PlateCarree)):
                    warnings.warn(f'Cannot add gridline labels on cartopy {self.projection} projection.')
                    labels = False
                # Grid label toggling
                if labels:
                    gl.ylabels_left   = latlabels[0]
                    gl.ylabels_right  = latlabels[1]
                    gl.xlabels_bottom = lonlabels[2]
                    gl.xlabels_top    = lonlabels[3]
                # Turn off gridlabels
                elif labels is not None:
                    gl.ylabels_left   = False
                    gl.ylabels_right  = False
                    gl.xlabels_bottom = False
                    gl.xlabels_top    = False
            # Turn off gridlines
            elif grid is not None:
                gl.xlines = False
                gl.ylines = False
                gl.ylabels_left   = False
                gl.ylabels_right  = False
                gl.xlabels_bottom = False
                gl.xlabels_top    = False

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
            reso = rc.get('reso')
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
            for name,args in features.items():
                # Get feature
                # TODO: Editing existing natural features? Creating natural
                # features at init time and hiding them?
                if not rc.get(name): # toggled
                    continue
                if getattr(self, f'_{name}', None): # already drawn
                    continue
                feat = cfeature.NaturalEarthFeature(*args, reso)
                # Customize
                # For 'lines', need to specify edgecolor and facecolor individually
                # See: https://github.com/SciTools/cartopy/issues/803
                kw = rc.category(name, cache=False)
                if name in ('coast', 'rivers', 'borders', 'innerborders'):
                    kw['edgecolor'] = kw.pop('color')
                    kw['facecolor'] = 'none'
                else:
                    kw['linewidth'] = 0
                if name in ('ocean',):
                    kw['zorder'] = 0.5 # below everything!
                self.add_feature(feat, **kw)
                setattr(self, '_' + name, feat)

            # Update patch
            kw_face = rc.fill({
                'facecolor': 'geoaxes.facecolor'
                })
            patch_kw = patch_kw or {}
            kw_face.update(patch_kw)
            self.background_patch.update(kw_face)
            kw_edge = rc.fill({
                'edgecolor': 'geoaxes.edgecolor',
                'linewidth': 'geoaxes.linewidth'
                })
            self.outline_patch.update(kw_edge)

            # Pass stuff to parent formatter, e.g. title and abc labeling
            super().format(**kwargs)

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Draw gridliner objects to tight bounding box algorithm will
        incorporate gridliner labels."""
        # Matplotlib approaches this problem in the same way, by duplicating
        # post-processing steps in both Axes.draw and Axes.get_tightbbox
        if self.get_autoscale_on() and self.ignore_existing_data_limits:
            self.autoscale_view()
        if self.background_patch.reclip:
            clipped_path = self.background_patch.orig_path.clip_to_bbox(
                self.viewLim)
            self.background_patch._path = clipped_path
        self.apply_aspect()
        for gl in self._gridliners:
            try: # new versions only
                gl._draw_gridliner(background_patch=self.background_patch,
                                renderer=renderer)
            except TypeError:
                gl._draw_gridliner(background_patch=self.background_patch)
        self._gridliners = []
        return super().get_tightbbox(renderer, *args, **kwargs)

class BasemapAxes(ProjectionAxes):
    """Axes subclass for plotting `~mpl_toolkits.basemap` projections. The
    `~mpl_toolkits.basemap.Basemap` projection instance is added as
    the `map_projection` attribute, but this is all abstracted away -- you can use
    `~matplotlib.axes.Axes` methods like `~matplotlib.axes.Axes.plot` and
    `~matplotlib.axes.Axes.contour` with your raw longitude-latitude data."""
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
            Passed to `Axes.__init__`.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.proj`
        """
        # Map boundary notes
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
        self.projection = map_projection
        self.boundary = None
        self._hasrecurred = False # use this so we can override plotting methods
        self._mapboundarydrawn = None
        self._latmax = None
        self._latlines = None
        self._lonlines = None
        self._lonlines_values = None
        self._latlines_values = None
        super().__init__(*args, **kwargs)

    def __getattribute__(self, attr, *args):
        """Applies the `~proplot.wrappers.cmap_wrapper`, `~proplot.wrappers.cycle_wrapper`,
        `~proplot.wrappers.enforce_centers`, `~proplot.wrappers.enforce_edges`,
        `~proplot.wrappers.basemap_gridfix`, `~proplot.wrappers.basemap_latlon`,
        `~proplot.wrappers.plot_wrapper`, `~proplot.wrappers.scatter_wrapper`,
        `~proplot.wrappers.fill_between_wrapper`, and `~proplot.wrappers.fill_betweenx_wrapper`
        wrappers. Also wraps all plotting methods with the hidden ``_basemap_call``
        and ``_no_recurse`` wrappers. Respectively, these call methods on
        the `~mpl_toolkits.basemap.Basemap` instance and prevent recursion
        issues arising from internal `~mpl_toolkits.basemap` calls to the
        axes methods."""
        # WARNING: Never ever try to just make blanket methods on the Basemap
        # instance accessible from axes instance! Can of worms and had bunch of
        # weird errors! Just pick the ones you think user will want to use.
        obj = super().__getattribute__(attr, *args)
        if attr in wrappers.LATLON_METHODS or attr in wrappers.EDGES_METHODS \
                or attr in wrappers.CENTERS_METHODS:
            # Step 6) Call identically named Basemap object method
            obj = wrappers._basemap_call(self, obj)
            # Step 5) Color usage wrappers
            if attr in wrappers.CMAP_METHODS:
                obj = wrappers._cmap_wrapper(self, obj)
            elif attr in wrappers.CYCLE_METHODS:
                obj = wrappers._cycle_wrapper(self, obj)
            # Step 4) Fix coordinate grid
            if attr in wrappers.EDGES_METHODS or attr in wrappers.CENTERS_METHODS:
                obj = wrappers._basemap_gridfix(self, obj)
            # Step 3) Utilities
            if attr in wrappers.EDGES_METHODS:
                obj = wrappers._enforce_edges(self, obj)
            elif attr in wrappers.CENTERS_METHODS:
                obj = wrappers._enforce_centers(self, obj)
            # Step 2) Better default keywords
            if attr in wrappers.LATLON_METHODS:
                obj = wrappers._basemap_latlon(self, obj)
            # Step 1) Parse args input
            if attr in wrappers.D2_METHODS:
                obj = wrappers._autoformat_2d(self, obj)
            elif attr in wrappers.D1_METHODS:
                obj = wrappers._autoformat_1d(self, obj)
            # Step 0) Special wrappers
            if attr == 'plot':
                obj = wrappers._plot_wrapper(self, obj)
            elif attr == 'scatter':
                obj = wrappers._scatter_wrapper(self, obj)
            elif attr == 'fill_between':
                obj = wrappers._fill_between_wrapper(self, obj)
            elif attr == 'fill_betweenx':
                obj = wrappers._fill_betweenx_wrapper(self, obj)
            # Recursion fix at top level
            obj = wrappers._no_recurse(self, obj)
        return obj

    def format(self, *, patch_kw=None, **kwargs):
        # Docstring added at bottom
        context, kwargs = self.context(**kwargs)
        with context:
            (grid, latmax, lonlim, latlim, boundinglat,
                lonlocator, latlocator, labels, lonlabels, latlabels,
                kwargs) = self._projection_format_kwargs(**kwargs)
            if (lonlim is not None or latlim is not None or
                boundinglat is not None):
                warnings.warn('Got lonlim={lonlim}, latlim={latlim}, boundinglat={boundinglat}, but you cannot "zoom into" a basemap projection after creating it. Pass a proj_kw dictionary in your call to subplots, with any of the following basemap keywords: boundinglat, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, llcrnrx, llcrnry, urcrnrx, urcrnry, width, or height.')
            # Map boundary
            # * First have to *manually replace* the old boundary by just
            #   deleting the original one
            # * If boundary is drawn successfully should be able to call
            #   self.projection._mapboundarydrawn.set_visible(False) and
            #   edges/fill color disappear
            # * For now will enforce that map plots *always* have background
            #   whereas axes plots can have transparent background
            kw_face = rc.fill({
                'facecolor': 'geoaxes.facecolor'
                })
            patch_kw = patch_kw or {}
            kw_face.update(patch_kw)
            kw_edge = rc.fill({
                'linewidth': 'geoaxes.linewidth',
                'edgecolor': 'geoaxes.edgecolor'
                })
            self.axesPatch = self.patch # bugfix or something
            if self.projection.projection in self._proj_non_rectangular:
                self.patch.set_alpha(0) # make patch invisible
                if not self.projection._mapboundarydrawn:
                    p = self.projection.drawmapboundary(ax=self) # set fill_color to 'none' to make transparent
                else:
                    p = self.projection._mapboundarydrawn
                p.update({**kw_face, **kw_edge})
                p.set_rasterized(False) # not sure about this; might be rasterized
                p.set_clip_on(False) # so edges of *line* denoting boundary aren't cut off
                self.boundary = p # not sure why this one
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
                # Change from left/right/bottom/top to left/right/top/bottom
                if labels:
                    lonlabels[2:] = lonlabels[2:][::-1]
                else:
                    lonlabels = 4*[0]
                if labels:
                    latlabels[2:] = latlabels[2:][::-1]
                else:
                    latlabels = 4*[0]
                # Turn off old ones
                # NOTE: Need to redraw both lon and lat lines if one changed,
                # because latmax affects their extent!
                if (latlocator is not None or lonlocator is not None
                    or latmax is not None or labels):
                    if self._lonlines:
                        for pi in self._lonlines.values():
                            for obj in [i for j in pi for i in j]: # magic
                                obj.set_visible(False)
                    if self._latlines:
                        for pi in self._latlines.values():
                            for obj in [i for j in pi for i in j]: # magic
                                obj.set_visible(False)
                    if latmax is None:
                        latmax = self._latmax
                    if latlocator is None:
                        latlocator = self._latlines_values
                    if lonlocator is None:
                        lonlocator = self._lonlines_values
                # Draw new ones
                if lonlocator is not None:
                    p = self.projection.drawmeridians(lonlocator,
                        latmax=latmax, labels=lonlabels, ax=self)
                    for pi in p.values():
                        for obj in [i for j in pi for i in j]: # magic
                            if isinstance(obj, mtext.Text):
                                obj.update(tkw)
                            else:
                                obj.update(lkw)
                    self._lonlines = p
                if latlocator is not None:
                    p = self.projection.drawparallels(latlocator,
                        latmax=latmax, labels=latlabels, ax=self)
                    for pi in p.values(): # returns dict, where each one is tuple
                        # Tried passing clip_on to the below, but it does nothing; must set
                        # for lines created after the fact
                        for obj in [i for j in pi for i in j]: # magic
                            if isinstance(obj, mtext.Text):
                                obj.update(tkw)
                            else:
                                obj.update(lkw)
                    self._latlines = p

            # Geography
            # TODO: Allow setting the zorder.
            # NOTE: Also notable are drawcounties, blumarble, drawlsmask,
            # shadedrelief, and etopo methods.
            features = {
                'land':         'fillcontinents',
                'coast':        'drawcoastlines',
                'rivers':       'drawrivers',
                'borders':      'drawcountries',
                'innerborders': 'drawstates',
                }
            for name, method in features.items():
                if not rc.get(name): # toggled
                    continue
                if getattr(self, f'_{name}', None): # already drawn
                    continue
                kw = rc.category(name, cache=False)
                feat = getattr(self.projection, method)(ax=self)
                if isinstance(feat, (list,tuple)): # can return single artist or list of artists
                    for obj in feat:
                        obj.update(kw)
                else:
                    feat.update(kw)
                setattr(self, '_' + name, feat)

            # Pass stuff to parent formatter, e.g. title and abc labeling
            Axes.format(self, **kwargs)

# Docstring setup
_projection_format_docstring = """
Calls `Axes.format` and `Axes.context`, formats the meridian
and parallel labels, longitude and latitude map limits, geographic
features, and more.

Parameters
----------
labels : bool, optional
    Whether to draw longitude and latitude labels. Defaults to
    ``rc['geogrid.labels']``.
lonlabels, latlabels
    Whether to label longitudes and latitudes, and on which sides
    of the map. There are four different options:

    1. Boolean ``True``. Indicates left side for latitudes,
       bottom for longitudes.
    2. A string, e.g. ``'lr'`` or ``'bt'``.
    3. A boolean ``(left,right)`` tuple for longitudes,
       ``(bottom,top)`` for latitudes.
    4. A boolean ``(left,right,bottom,top)`` tuple as in the
       `~mpl_toolkits.basemap.Basemap.drawmeridians` and
       `~mpl_toolkits.basemap.Basemap.drawparallels` methods.

lonlim, latlim : (float, float), optional
    Longitude and latitude limits of projection, applied
    with `~cartopy.mpl.geoaxes.GeoAxes.set_extent`. For cartopy axes only.
boundinglat : float, optional
    The edge latitude for the circle bounding North Pole and
    South Pole-centered projections. For cartopy axes only.
grid : bool, optional
    Whether to add meridian and parallel gridlines.
latmax : float, optional
    The maximum absolute latitude for meridian gridlines. Defaults
    to ``rc['geogrid.latmax']``.
lonlines, latlines : float or list of float, optional
    If float, indicates the *spacing* of meridian and parallel gridlines.
    Otherwise, must be a list of floats indicating specific meridian and
    parallel gridlines to draw.
lonlocator, latlocator : optional
    Aliases for `lonlines`, `latlines`.
patch_kw : dict-like, optional
    Keyword arguments used to update the background patch object. You
    can use this, for example, to set background hatching with
    ``patch_kw={'hatch':'xxx'}``.
**kwargs
    Passed to `Axes.format` and `Axes.context`.
"""
CartopyAxes.format.__doc__ = _projection_format_docstring
BasemapAxes.format.__doc__ = _projection_format_docstring

# Register the projections
# TODO: Remove BasemapAxes!!! Cartopy will support gridline labels soon.
mproj.register_projection(PolarAxes)
mproj.register_projection(CartesianAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)

