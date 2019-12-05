#!/usr/bin/env python3
"""
The axes classes used for all ProPlot figures.
"""
import numpy as np
import warnings
import functools
from numbers import Integral
import matplotlib.projections as mproj
import matplotlib.axes as maxes
import matplotlib.dates as mdates
import matplotlib.scale as mscale
import matplotlib.text as mtext
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections
from . import utils, projs, axistools
from .utils import _notNone, units
from .rctools import rc, RC_NODOTSNAMES
from .wrappers import (
    _get_transform, _norecurse, _redirect,
    _add_errorbars, _bar_wrapper, _barh_wrapper, _boxplot_wrapper,
    _default_crs, _default_latlon, _default_transform, _cmap_changer,
    _cycle_changer, _fill_between_wrapper, _fill_betweenx_wrapper,
    _hist_wrapper, _plot_wrapper, _scatter_wrapper,
    _standardize_1d, _standardize_2d,
    _text_wrapper, _violinplot_wrapper,
    colorbar_wrapper, legend_wrapper,
)
try:
    from cartopy.mpl.geoaxes import GeoAxes
except ModuleNotFoundError:
    GeoAxes = object

__all__ = [
    'Axes',
    'BasemapAxes',
    'GeoAxes',
    'PolarAxes', 'ProjAxes',
    'XYAxes',
]

# Translator for inset colorbars and legends
ABC_STRING = 'abcdefghijklmnopqrstuvwxyz'
LOC_TRANSLATE = {
    None: None,
    'inset': 'best',
    'i': 'best',
    0: 'best',
    1: 'upper right',
    2: 'upper left',
    3: 'lower left',
    4: 'lower right',
    5: 'center left',
    6: 'center right',
    7: 'lower center',
    8: 'upper center',
    9: 'center',
    'l': 'left',
    'r': 'right',
    'b': 'bottom',
    't': 'top',
    'c': 'center',
    'ur': 'upper right',
    'ul': 'upper left',
    'll': 'lower left',
    'lr': 'lower right',
    'cr': 'center right',
    'cl': 'center left',
    'uc': 'upper center',
    'lc': 'lower center',
}
SIDE_TRANSLATE = {
    'l': 'left',
    'r': 'right',
    'b': 'bottom',
    't': 'top',
}


def _abc(i):
    """Function for a-b-c labeling, returns a...z...aa...zz...aaa...zzz."""
    if i < 26:
        return ABC_STRING[i]
    else:
        return _abc(i - 26) + ABC_STRING[i % 26]  # sexy sexy recursion


def _disable_decorator(msg):
    """
    Generate decorators that disable methods. Also sets __doc__ to None so
    that ProPlot fork of automodapi doesn't add these methods to the website
    documentation. Users can still call help(ax.method) because python looks
    for superclass method docstrings if a docstring is empty.
    """
    def decorator(func):
        @functools.wraps(func)
        def _wrapper(self, *args, **kwargs):
            raise RuntimeError(msg.format(func.__name__))
        _wrapper.__doc__ = None
        return _wrapper
    return decorator


class Axes(maxes.Axes):
    """Lowest-level axes subclass. Handles titles and axis
    sharing. Adds several new methods and overrides existing ones."""

    def __init__(self, *args, number=None,
                 sharex=0, sharey=0,
                 spanx=None, spany=None, alignx=None, aligny=None,
                 main=False,
                 **kwargs):
        """
        Parameters
        ----------
        number : int
            The subplot number, used for a-b-c labelling (see
            `~Axes.format`).
        sharex, sharey : {3, 2, 1, 0}, optional
            The "axis sharing level" for the *x* axis, *y* axis, or both
            axes. See `~proplot.subplots.subplots` for details.
        spanx, spany : bool, optional
            Boolean toggle for whether spanning labels are enabled for the
            *x* and *y* axes. See `~proplot.subplots.subplots` for details.
        alignx, aligny : bool, optional
            Boolean toggle for whether aligned axis labels are enabled for the
            *x* and *y* axes. See `~proplot.subplots.subplots` for details.
        main : bool, optional
            Used internally, indicates whether this is a "main axes" rather
            than a twin, panel, or inset axes.

        See also
        --------
        :py:obj:`matplotlib.axes.Axes`,
        :py:obj:`XYAxes`,
        :py:obj:`PolarAxes`,
        :py:obj:`ProjAxes`
        """
        # Call parent
        super().__init__(*args, **kwargs)
        # Properties
        self._number = number  # for abc numbering
        self._abc_loc = None
        self._abc_text = None
        self._titles_dict = {}  # title text objects and their locations
        self._title_loc = None  # location of main title
        # so we can copy to top panel
        self._title_pad = rc.get('axes.titlepad')
        self._title_above_panel = True  # TODO: add rc prop?
        # Children and related properties
        self._bpanels = []
        self._tpanels = []
        self._lpanels = []
        self._rpanels = []
        self._tightbbox = None  # bounding boxes are saved
        self._panel_side = None
        self._panel_share = False  # True when "filled" with cbar/legend
        self._panel_parent = None
        self._panel_filled = False  # True when "filled" with cbar/legend
        self._inset_parent = None
        self._inset_zoom = False
        self._inset_zoom_data = None
        self._alty_child = None
        self._altx_child = None
        self._alty_parent = None
        self._altx_parent = None
        self._auto_colorbar = {}  # stores handles and kwargs for auto colorbar
        self._auto_legend = {}
        # Text labels
        # TODO: Add text labels as panels instead of as axes children?
        coltransform = mtransforms.blended_transform_factory(
            self.transAxes, self.figure.transFigure)
        rowtransform = mtransforms.blended_transform_factory(
            self.figure.transFigure, self.transAxes)
        self._llabel = self.text(
            0.05, 0.5, '', va='center', ha='right', transform=rowtransform)
        self._rlabel = self.text(
            0.95, 0.5, '', va='center', ha='left', transform=rowtransform)
        self._blabel = self.text(
            0.5, 0.05, '', va='top', ha='center', transform=coltransform)
        self._tlabel = self.text(
            0.5, 0.95, '', va='bottom', ha='center', transform=coltransform)
        # Shared and spanning axes
        if main:
            self.figure._axes_main.append(self)
        self._spanx_on = spanx
        self._spany_on = spany
        self._alignx_on = alignx
        self._aligny_on = aligny
        self._sharex_level = sharex
        self._sharey_level = sharey
        self.format(mode=1)  # mode == 1 applies the rcShortParams

    def _draw_auto_legends_colorbars(self):
        """Generate automatic legends and colorbars. Wrapper funcs
        let user add handles to location lists with successive calls to
        make successive calls to plotting commands."""
        for loc, (handles, kwargs) in self._auto_colorbar.items():
            self.colorbar(handles, **kwargs)
        for loc, (handles, kwargs) in self._auto_legend.items():
            self.legend(handles, **kwargs)
        self._auto_legend = {}
        self._auto_colorbar = {}

    def _get_side_axes(self, side):
        """Returns axes whose left, right, top, or bottom side abutts
        against the same row or column as this axes."""
        s = side[0]
        if s not in 'lrbt':
            raise ValueError(f'Invalid side {side!r}.')
        if not hasattr(self, 'get_subplotspec'):
            return [self]
        x = ('x' if s in 'lr' else 'y')
        idx = (0 if s in 'lt' else 1)  # which side of range to test
        coord = self._range_gridspec(x)[idx]  # side for a particular axes
        axs = [ax for ax in self.figure._axes_main
               if ax._range_gridspec(x)[idx] == coord]
        if not axs:
            return [self]
        else:
            return axs

    def _get_extent_axes(self, x):
        """Returns axes whose horizontal or vertical extent in the main
        gridspec matches the horizontal or vertical extend of this axes.
        Also sorts the list so the leftmost or bottommost axes is at the
        start of the list."""
        if not hasattr(self, 'get_subplotspec'):
            return [self]
        y = ('y' if x == 'x' else 'x')
        idx = (0 if x == 'x' else 1)
        argfunc = (np.argmax if x == 'x' else np.argmin)
        irange = self._range_gridspec(x)
        axs = [ax for ax in self.figure._axes_main
               if ax._range_gridspec(x) == irange]
        if not axs:
            return [self]
        else:
            pax = axs.pop(argfunc([ax._range_gridspec(y)[idx] for ax in axs]))
            return [pax, *axs]

    def _get_title_props(self, abc=False, loc=None):
        """Returns standardized location name, position keyword arguments, and
        setting keyword arguments for the relevant title or a-b-c label at
        location `loc`."""
        # Props
        # NOTE: Sometimes we load all properties from rc object, sometimes
        # just changed ones. This is important if e.g. user calls in two
        # lines ax.format(titleweight='bold') then ax.format(title='text'),
        # don't want to override custom setting with rc default setting.
        def props(cache):
            return rc.fill({
                'fontsize': f'{prefix}.size',
                'weight': f'{prefix}.weight',
                'color': f'{prefix}.color',
                'border': f'{prefix}.border',
                'linewidth': f'{prefix}.linewidth',
                'fontfamily': 'font.family',
            }, cache=cache)

        # Location string and position coordinates
        cache = True
        prefix = 'abc' if abc else 'title'
        loc = _notNone(loc, rc[f'{prefix}.loc'])
        iloc = getattr(self, '_' + ('abc' if abc else 'title') + '_loc')  # old
        if loc is None:
            loc = iloc
        elif iloc is not None and loc != iloc:
            cache = False

        # Above axes
        loc = LOC_TRANSLATE.get(loc, loc)
        if loc in ('top', 'bottom'):
            raise ValueError(f'Invalid title location {loc!r}.')
        elif loc in ('left', 'right', 'center'):
            kw = props(cache)
            kw.pop('border', None)  # no border for titles outside axes
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
            if loc in ('upper center', 'lower center'):
                x, ha = 0.5, 'center'
            elif loc in ('upper left', 'lower left'):
                xpad = rc.get('axes.titlepad') / (72 * width)
                x, ha = 1.5 * xpad, 'left'
            elif loc in ('upper right', 'lower right'):
                xpad = rc.get('axes.titlepad') / (72 * width)
                x, ha = 1 - 1.5 * xpad, 'right'
            else:
                raise ValueError(f'Invalid title or abc "loc" {loc}.')
            if loc in ('upper left', 'upper right', 'upper center'):
                ypad = rc.get('axes.titlepad') / (72 * height)
                y, va = 1 - 1.5 * ypad, 'top'
            elif loc in ('lower left', 'lower right', 'lower center'):
                ypad = rc.get('axes.titlepad') / (72 * height)
                y, va = 1.5 * ypad, 'bottom'
            obj = self.text(x, y, '', ha=ha, va=va, transform=self.transAxes)
            obj.set_transform(self.transAxes)
        return loc, obj, kw

    @staticmethod
    def _loc_translate(loc, **kwargs):
        """Translates location string `loc` into a standardized form."""
        if loc is True:
            loc = 'r'  # for on-the-fly colorbars and legends
        elif isinstance(loc, (str, Integral)):
            loc = LOC_TRANSLATE.get(loc, loc)
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
        """Gets the column or row range for the axes."""
        subplotspec = self.get_subplotspec()
        if x == 'x':
            _, _, _, _, col1, col2 = subplotspec.get_active_rows_columns()
            return col1, col2
        else:
            _, _, row1, row2, _, _ = subplotspec.get_active_rows_columns()
            return row1, row2

    def _range_tightbbox(self, x):
        """Gets span of tight bounding box, including twin axes and panels
        which are not considered real children and so aren't ordinarily
        included in the tight bounding box calc.
        `~proplot.axes.Axes.get_tightbbox` caches tight bounding boxes when
        `~Figure.get_tightbbox` is called."""
        # TODO: Better resting for axes visibility
        bbox = self._tightbbox
        if bbox is None:
            return np.nan, np.nan
        if x == 'x':
            return bbox.xmin, bbox.xmax
        else:
            return bbox.ymin, bbox.ymax

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
        for key in ('color', 'fontproperties'):  # TODO: add to this?
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
        if self._panel_side == 'top' and self._panel_parent:
            ax, taxs = self._panel_parent, [self]
        else:
            ax, taxs = self, self._tpanels
        if not taxs or not ax._title_above_panel:
            tax = ax
        else:
            tax = taxs[0]
            tax._title_pad = ax._title_pad
            for loc, obj in ax._titles_dict.items():
                if not obj.get_text() or loc not in (
                        'left', 'center', 'right'):
                    continue
                kw = {}
                loc, tobj, _ = tax._get_title_props(loc=loc)
                for key in ('text', 'color', 'fontproperties'):  # add to this?
                    kw[key] = getattr(obj, 'get_' + key)()
                tobj.update(kw)
                tax._titles_dict[loc] = tobj
                obj.set_text('')

        # Push title above tick marks -- this is known matplotlib problem,
        # but especially annoying with top panels!
        # TODO: Make sure this is robust. Seems 'default' is returned usually
        # when tick label sides is actually *both*. Also makes sure axis is
        # visible; if not, this is a filled cbar/legend, no padding needed
        pad = 0
        pos = tax.xaxis.get_ticks_position()
        labs = tax.xaxis.get_ticklabels()
        if pos == 'default' or (pos == 'top' and not len(labs)) or (
                pos == 'unknown' and tax._panel_side == 'top'
                and not len(labs) and tax.xaxis.get_visible()):
            pad = tax.xaxis.get_tick_padding()
        tax._set_title_offset_trans(self._title_pad + pad)

    def _sharex_setup(self, sharex, level):
        """Sets up panel axis sharing."""
        if level not in range(4):
            raise ValueError(
                'Level can be 0 (share nothing), '
                '1 (do not share limits, just hide axis labels), '
                '2 (share limits, but do not hide tick labels), or '
                '3 (share limits and hide tick labels). Got {level}.')
        # enforce, e.g. if doing panel sharing
        self._sharex_level = max(self._sharex_level, level)
        self._share_short_axis(sharex, 'l', level)
        self._share_short_axis(sharex, 'r', level)
        self._share_long_axis(sharex, 'b', level)
        self._share_long_axis(sharex, 't', level)

    def _sharey_setup(self, sharey, level):
        """Sets up panel axis sharing."""
        if level not in range(4):
            raise ValueError(
                'Level can be 0 (share nothing), '
                '1 (do not share limits, just hide axis labels), '
                '2 (share limits, but do not hide tick labels), or '
                '3 (share limits and hide tick labels). Got {level}.')
        self._sharey_level = max(self._sharey_level, level)
        self._share_short_axis(sharey, 'b', level)
        self._share_short_axis(sharey, 't', level)
        self._share_long_axis(sharey, 'l', level)
        self._share_long_axis(sharey, 'r', level)

    def _share_setup(self):
        """Applies axis sharing for axes that share the same horizontal or
        vertical extent, and for their panels."""
        # Panel axes sharing, between main subplot and its panels
        # Top and bottom
        def shared(paxs):
            return [
                pax for pax in paxs if not pax._panel_filled
                and pax._panel_share]
        if not self._panel_side:  # this is a main axes
            bottom = self
            paxs = shared(self._bpanels)
            if paxs:
                bottom = paxs[-1]
                for iax in (self, *paxs[:-1]):
                    # parent is *bottom-most* panel
                    iax._sharex_setup(bottom, 3)
            paxs = shared(self._tpanels)
            for iax in paxs:
                iax._sharex_setup(bottom, 3)
            # Left and right
            left = self
            paxs = shared(self._lpanels)
            if paxs:
                left = paxs[0]
                for iax in (*paxs[1:], self):
                    iax._sharey_setup(left, 3)  # parent is *bottom-most* panel
            paxs = shared(self._rpanels)
            for iax in paxs:
                iax._sharey_setup(left, 3)

        # Main axes, sometimes overrides panel axes sharing
        # TODO: This can get very repetitive, but probably minimal impact
        # on performance?
        # Share x axes
        parent, *children = self._get_extent_axes('x')
        for child in children:
            child._sharex_setup(parent, parent._sharex_level)
        # Share y axes
        parent, *children = self._get_extent_axes('y')
        for child in children:
            child._sharey_setup(parent, parent._sharey_level)

    def _share_short_axis(self, share, side, level):
        """Share the "short" axes of panels along a main subplot with panels
        along an external subplot."""
        if share is None or self._panel_side:  # not None
            return
        s = side[0]
        axis = 'x' if s in 'lr' else 'y'
        caxs = getattr(self, '_' + s + 'panels')
        paxs = getattr(share, '_' + s + 'panels')
        caxs = [pax for pax in caxs if not pax._panel_filled]
        paxs = [pax for pax in paxs if not pax._panel_filled]
        for cax, pax in zip(caxs, paxs):  # may be uneven
            getattr(cax, '_share' + axis + '_setup')(pax, level)

    def _share_long_axis(self, share, side, level):
        """Share the "long" axes of panels along a main subplot with the
        axis from an external subplot."""
        # NOTE: We do not check _panel_share because that only controls
        # sharing with main subplot, not other subplots
        if share is None or self._panel_side:
            return
        s = side[0]
        axis = 'x' if s in 'tb' else 'y'
        paxs = getattr(self, '_' + s + 'panels')
        paxs = [pax for pax in paxs if not pax._panel_filled]
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share, level)

    def _update_axislabels(self, x='x', **kwargs):
        """Apply axis labels to the relevant shared axis. If spanning
        labels are toggled this keeps the labels synced for all subplots in
        the same row or column. Label positions will be adjusted at draw-time
        with figure._align_axislabels."""
        if x not in 'xy':
            return
        # Update label on this axes
        axis = getattr(self, x + 'axis')
        axis.label.update(kwargs)
        kwargs.pop('color', None)

        # Defer to parent (main) axes if possible, then get the axes
        # shared by that parent
        ax = self._panel_parent or self
        ax = getattr(ax, '_share' + x) or ax

        # Apply to spanning axes and their panels
        axs = [ax]
        if getattr(ax, '_span' + x + '_on'):
            s = axis.get_label_position()[0]
            if s in 'lb':
                axs = ax._get_side_axes(s)
        for ax in axs:
            getattr(ax, x + 'axis').label.update(kwargs)  # apply to main axes
            pax = getattr(ax, '_share' + x)
            if pax is not None:  # apply to panel?
                getattr(pax, x + 'axis').label.update(kwargs)

    def _update_title(self, obj, **kwargs):
        """Redraws title if updating with input keyword args failed."""
        # Try to just return updated object, redraw may be necessary
        # WARNING: Making text instances invisible seems to mess up tight
        # bounding box calculations and cause other issues. Just reset text.
        keys = ('border', 'lw', 'linewidth', 'bordercolor', 'invert')
        kwextra = {key: value for key, value in kwargs.items() if key in keys}
        kwargs = {key: value for key,
                  value in kwargs.items() if key not in keys}
        obj.update(kwargs)
        if kwextra:
            obj.set_text('')
        else:
            return obj
        # Get properties from old object
        for key in ('ha', 'va', 'color', 'transform', 'fontproperties'):
            kwextra[key] = getattr(obj, 'get_' + key)()  # copy over attrs
        text = kwargs.pop('text', obj.get_text())
        x, y = kwargs.pop('position', (None, None))
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
            `~proplot.rctools.rc` object.
        **kwargs
            Any of three options:

            * A keyword arg for `Axes.format`, `XYAxes.format`,
              or `ProjAxes.format`.
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
            `~proplot.rctools.rc` object will ignore
            `rcParams <https://matplotlib.org/users/customizing.html>`__.
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
        kw = {}  # for format
        rc_kw = rc_kw or {}
        for key, value in kwargs.items():
            key_fixed = RC_NODOTSNAMES.get(key, None)
            if key_fixed is None:
                kw[key] = value
            else:
                rc_kw[key_fixed] = value
        rc._getitem_mode = 0  # might still be non-zero if had error
        # Return "context object", which is just the configurator itself
        # primed for use in a "with" statement
        return rc.context(rc_kw, mode=mode), kw

    def format(
            self, *, title=None, top=None,
            figtitle=None, suptitle=None, rowlabels=None, collabels=None,
            leftlabels=None, rightlabels=None,
            toplabels=None, bottomlabels=None,
            llabels=None, rlabels=None, tlabels=None, blabels=None,
            **kwargs):
        """
        Called by `XYAxes.format`, `ProjAxes.format`, and
        `PolarAxes.format`. Formats the axes title(s), the a-b-c label, row
        and column labels, and the figure title.

        Parameters
        ----------
        title : str, optional
            The axes title.
        abc : bool, optional
            Whether to apply "a-b-c" subplot labelling based on the
            ``number`` attribute. If ``number`` is >26, the labels will loop
            around to a, ..., z, aa, ..., zz, aaa, ..., zzz, ... Default is
            :rc:`abc`.
        abcstyle : str, optional
            String denoting the format of a-b-c labels containing the character
            ``a`` or ``A``. ``'a'`` is the default, but e.g. ``'a.'``,
            ``'a)'``, or ``'A'`` might also be desirable. Default is
            :rc:`abc.style`.
        abcloc, titleloc : str, optional
            Strings indicating the location for the a-b-c label and
            main title. The following locations keys are valid. Defaults are
            :rc:`abc.loc` and :rc:`title.loc`.

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
            of artists plotted inside the axes. Defaults are
            :rc:`abc.border` and :rc:`title.border`
        ltitle, rtitle, ultitle, uctitle, urtitle, lltitle, lctitle, lrtitle \
: str, optional
            Axes titles in particular positions. This lets you specify multiple
            "titles" for each subplots. See the `abcloc` keyword.
        top : bool, optional
            Whether to try to put title and a-b-c label above the top subplot
            panel (if it exists), or to always put them on the main subplot.
            Default is ``True``.
        rowlabels, colllabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`.
        llabels, tlabels, rlabels, blabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`, `rightlabels`,
            `bottomlabels`.
        leftlabels, toplabels, rightlabels, bottomlabels : list of str, \
optional
            The subplot row and column labels. If list, length must match
            the number of subplots on the left, top, right, or bottom edges
            of the figure.
        figtitle, suptitle : str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.
            This is an improvement on matplotlib's "super" title, which just
            centers the text between figure edges.

        Note
        ----
        The `abc`, `abcstyle`, `abcloc`, and `titleloc` keyword arguments
        are actually rc configuration settings that are temporarily
        changed by the call to `~Axes.context`. They are documented here
        because it is extremely common to change them with `~Axes.format`.
        They also appear in the tables in the `~proplot.rctools` documention.

        See also
        --------
        :py:obj:`Axes.context`,
        :py:obj:`XYAxes.format`,
        :py:obj:`ProjAxes.format`,
        :py:obj:`PolarAxes.format`,
        """
        # Figure patch (for some reason needs to be re-asserted even if
        # declared before figure is drawn)
        kw = rc.fill({'facecolor': 'figure.facecolor'})
        self.figure.patch.update(kw)
        if top is not None:
            self._title_above_panel = top
        pad = rc['axes.titlepad']
        if pad is not None:
            self._set_title_offset_trans(pad)
            self._title_pad = pad

        # Super title
        # NOTE: These are actually *figure-wide* settings, but that line
        # gets blurred where we have shared axes, spanning labels, and
        # whatnot. May result in redundant assignments if formatting more than
        # one axes, but operations are fast so some redundancy is nbd.
        # NOTE: Below workaround prevents changed *figure-wide* settings
        # from getting overwritten when user makes a new axes.
        fig = self.figure
        suptitle = _notNone(figtitle, suptitle, None,
                            names=('figtitle', 'suptitle'))
        if len(fig._axes_main) > 1 and rc._getitem_mode == 1:
            kw = {}
        else:
            kw = rc.fill({
                'fontsize': 'suptitle.size',
                'weight': 'suptitle.weight',
                'color': 'suptitle.color',
                'fontfamily': 'font.family'
            })
        if suptitle or kw:
            fig._update_figtitle(suptitle, **kw)
        # Labels
        llabels = _notNone(rowlabels, leftlabels, llabels,
                           None, names=('rowlabels', 'leftlabels', 'llabels'))
        tlabels = _notNone(collabels, toplabels, tlabels,
                           None, names=('collabels', 'toplabels', 'tlabels'))
        rlabels = _notNone(rightlabels, rlabels, None,
                           names=('rightlabels', 'rlabels'))
        blabels = _notNone(bottomlabels, blabels, None,
                           names=('bottomlabels', 'blabels'))
        for side, labels in zip(
                ('left', 'right', 'top', 'bottom'),
                (llabels, rlabels, tlabels, blabels)):
            kw = rc.fill({
                'fontsize': side + 'label.size',
                'weight': side + 'label.weight',
                'color': side + 'label.color',
                'fontfamily': 'font.family'
            })
            if labels or kw:
                fig._update_labels(self, side, labels, **kw)

        # A-b-c labels
        titles_dict = self._titles_dict
        if not self._panel_side:
            # Location and text
            abcstyle = rc['abc.style']  # changed or running format first time?
            if 'abcformat' in kwargs:  # super sophisticated deprecation system
                abcstyle = kwargs.pop('abcformat')
                warnings.warn(
                    f'rc setting "abcformat" is deprecated. '
                    f'Please use "abcstyle".')
            if abcstyle and self.number is not None:
                if not isinstance(abcstyle, str) or (
                        abcstyle.count('a') != 1 and abcstyle.count('A') != 1):
                    raise ValueError(
                        f'Invalid abcstyle {abcstyle!r}. '
                        'Must include letter "a" or "A".')
                abcedges = abcstyle.split('a' if 'a' in abcstyle else 'A')
                text = abcedges[0] + _abc(self.number - 1) + abcedges[-1]
                if 'A' in abcstyle:
                    text = text.upper()
                self._abc_text = text
            # Apply new settings
            # Also if a-b-c label was moved, remove previous one and update
            # text on new one, in case self._abc_text has not changed.
            loc, obj, kw = self._get_title_props(abc=True)
            iloc = self._abc_loc
            obj = self._update_title(obj, **kw)
            titles_dict[loc] = obj
            if iloc is not None and loc != iloc:
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
            for iloc, iobj in titles_dict.items():
                if iloc is self._abc_loc:
                    continue
                titles_dict[iloc] = self._update_title(iobj, **kw)
        # Workflow 2, want this to come first so workflow 1 gets priority
        for ikey, ititle in kwargs.items():
            if not ikey[-5:] == 'title':
                raise TypeError(
                    f'format() got an unexpected keyword argument {ikey!r}.')
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
        self._title_loc = loc  # assigns default loc on first run
        if title is not None:
            kw['text'] = title
        if kw:
            titles_dict[loc] = self._update_title(obj, **kw)

    def area(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.fill_between`."""
        # NOTE: *Cannot* assign area = axes.Axes.fill_between because the
        # wrapper won't be applied and for some reason it messes up
        # automodsumm, which tries to put the matplotlib docstring on website
        return self.fill_between(*args, **kwargs)

    def areax(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.fill_betweenx`."""
        return self.fill_betweenx(*args, **kwargs)

    def boxes(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.boxplot`."""
        return self.boxplot(*args, **kwargs)

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
            The colorbar location. Default is :rc:`colorbar.loc`. The
            following location keys are valid.

            ==================  ==================================
            Location            Valid keys
            ==================  ==================================
            outer left          ``'left'``, ``'l'``
            outer right         ``'right'``, ``'r'``
            outer bottom        ``'bottom'``, ``'b'``
            outer top           ``'top'``, ``'t'``
            default inset       ``'inset'``, ``'i'``, ``0``
            upper right inset   ``'upper right'``, ``'ur'``, ``1``
            upper left inset    ``'upper left'``, ``'ul'``, ``2``
            lower left inset    ``'lower left'``, ``'ll'``, ``3``
            lower right inset   ``'lower right'``, ``'lr'``, ``4``
            ==================  ==================================

        pad : float or str, optional
            The space between the axes edge and the colorbar. For inset
            colorbars only. Units are interpreted by `~proplot.utils.units`.
            Default is :rc:`colorbar.axespad`.
        length : float or str, optional
            The colorbar length. For outer colorbars, units are relative to the
            axes width or height. Default is :rc:`colorbar.length`. For inset
            colorbars, units are interpreted by `~proplot.utils.units`. Default
            is :rc:`colorbar.insetlength`.
        width : float or str, optional
            The colorbar width. Units are interpreted by
            `~proplot.utils.units`. Default is :rc:`colorbar.width` or
            :rc:`colorbar.insetwidth`.
        space : float or str, optional
            The space between the colorbar and the main axes. For outer
            colorbars only. Units are interpreted by `~proplot.utils.units`.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, defaut is :rc:`subplots.panelspace`.
        frame, frameon : bool, optional
            Whether to draw a frame around inset colorbars, just like
            `~matplotlib.axes.Axes.legend`.
            Default is :rc:`colorbar.frameon`.
        alpha, linewidth, edgecolor, facecolor : optional
            Transparency, edge width, edge color, and face color for the frame
            around the inset colorbar. Default is
            :rc:`colorbar.framealpha`, :rc:`axes.linewidth`,
            :rc:`axes.edgecolor`, and :rc:`axes.facecolor`,
            respectively.
        **kwargs
            Passed to `~proplot.wrappers.colorbar_wrapper`.
        """
        # TODO: add option to pad inset away from axes edge!
        kwargs.update({'edgecolor': edgecolor, 'linewidth': linewidth})
        loc = _notNone(loc, rc['colorbar.loc'])
        loc = self._loc_translate(loc)
        if loc == 'best':  # a white lie
            loc = 'lower right'
        if not isinstance(loc, str):  # e.g. 2-tuple or ndarray
            raise ValueError(f'Invalid colorbar location {loc!r}.')

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            return ax.colorbar(loc='_fill', *args, length=length, **kwargs)

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

            # Draw colorbar with arbitrary length relative to full length
            # of panel
            side = self._panel_side
            length = _notNone(length, rc['colorbar.length'])
            subplotspec = self.get_subplotspec()
            if length <= 0 or length > 1:
                raise ValueError(
                    f'Panel colorbar length must satisfy 0 < length <= 1, '
                    f'got length={length!r}.')
            if side in ('bottom', 'top'):
                gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows=1, ncols=3, wspace=0,
                    subplot_spec=subplotspec,
                    width_ratios=((1 - length) / 2, length, (1 - length) / 2),
                )
                subplotspec = gridspec[1]
            else:
                gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows=3, ncols=1, hspace=0,
                    subplot_spec=subplotspec,
                    height_ratios=((1 - length) / 2, length, (1 - length) / 2),
                )
                subplotspec = gridspec[1]
            with self.figure._unlock():
                ax = self.figure.add_subplot(subplotspec, projection=None)
            if ax is self:
                raise ValueError(f'Uh oh.')
            self.add_child_axes(ax)

            # Location
            if side in ('bottom', 'top'):
                outside, inside = 'bottom', 'top'
                if side == 'top':
                    outside, inside = inside, outside
                ticklocation = outside
                orientation = 'horizontal'
            else:
                outside, inside = 'left', 'right'
                if side == 'right':
                    outside, inside = inside, outside
                ticklocation = outside
                orientation = 'vertical'

            # Keyword args and add as child axes
            orient = kwargs.get('orientation', None)
            if orient is not None and orient != orientation:
                warnings.warn(f'Overriding input orientation={orient!r}.')
            ticklocation = kwargs.pop('tickloc', None) or ticklocation
            ticklocation = kwargs.pop('ticklocation', None) or ticklocation
            kwargs.update({'orientation': orientation,
                           'ticklocation': ticklocation})

        # Inset colorbar
        else:
            # Default props
            cbwidth, cblength = width, length
            width, height = self.get_size_inches()
            extend = units(_notNone(
                kwargs.get('extendsize', None), rc['colorbar.insetextend']))
            cbwidth = units(_notNone(
                cbwidth, rc['colorbar.insetwidth'])) / height
            cblength = units(_notNone(
                cblength, rc['colorbar.insetlength'])) / width
            pad = units(_notNone(pad, rc['colorbar.axespad']))
            xpad, ypad = pad / width, pad / height

            # Get location in axes-relative coordinates
            # Bounds are x0, y0, width, height in axes-relative coordinates
            xspace = rc['xtick.major.size'] / 72
            if kwargs.get('label', ''):
                xspace += 2.4 * rc['font.size'] / 72
            else:
                xspace += 1.2 * rc['font.size'] / 72
            xspace /= height  # space for labels
            if loc == 'upper right':
                bounds = (1 - xpad - cblength, 1 - ypad - cbwidth)
                fbounds = (1 - 2 * xpad - cblength,
                           1 - 2 * ypad - cbwidth - xspace)
            elif loc == 'upper left':
                bounds = (xpad, 1 - ypad - cbwidth)
                fbounds = (0, 1 - 2 * ypad - cbwidth - xspace)
            elif loc == 'lower left':
                bounds = (xpad, ypad + xspace)
                fbounds = (0, 0)
            elif loc == 'lower right':
                bounds = (1 - xpad - cblength, ypad + xspace)
                fbounds = (1 - 2 * xpad - cblength, 0)
            else:
                raise ValueError(f'Invalid colorbar location {loc!r}.')
            bounds = (bounds[0], bounds[1], cblength, cbwidth)
            fbounds = (fbounds[0], fbounds[1],
                       2 * xpad + cblength, 2 * ypad + cbwidth + xspace)

            # Make frame
            # NOTE: We do not allow shadow effects or fancy edges effect.
            # Also keep zorder same as with legend.
            frameon = _notNone(
                frame, frameon, rc['colorbar.frameon'],
                names=('frame', 'frameon'))
            if frameon:
                # Make patch object
                xmin, ymin, width, height = fbounds
                patch = mpatches.Rectangle(
                    (xmin, ymin), width, height,
                    snap=True, zorder=4, transform=self.transAxes)
                # Update patch props
                alpha = _notNone(alpha, rc['colorbar.framealpha'])
                linewidth = _notNone(linewidth, rc['axes.linewidth'])
                edgecolor = _notNone(edgecolor, rc['axes.edgecolor'])
                facecolor = _notNone(facecolor, rc['axes.facecolor'])
                patch.update({'alpha': alpha, 'linewidth': linewidth,
                              'edgecolor': edgecolor, 'facecolor': facecolor})
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
                warnings.warn(
                    f'Orientation for inset colorbars must be horizontal, '
                    f'ignoring orient={orient!r}.')
            ticklocation = kwargs.pop('tickloc', None)
            ticklocation = kwargs.pop('ticklocation', None) or ticklocation
            if ticklocation is not None and ticklocation != 'bottom':
                warnings.warn(
                    f'Inset colorbars can only have ticks on the bottom.')
            kwargs.update({'orientation': 'horizontal',
                           'ticklocation': 'bottom'})
            kwargs.setdefault('maxn', 5)
            kwargs.setdefault('extendsize', extend)

        # Generate colorbar
        return colorbar_wrapper(ax, *args, **kwargs)

    def legend(self, *args, loc=None, width=None, space=None, **kwargs):
        """
        Adds an *inset* legend or *outer* legend along the edge of the axes.
        See `~proplot.wrappers.legend_wrapper` for details.

        Parameters
        ----------
        loc : int or str, optional
            The legend location or panel location. The following location keys
            are valid. Note that if a panel does not exist, it will be
            generated on-the-fly.

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            left panel          ``'left'``, ``'l'``
            right panel         ``'right'``, ``'r'``
            bottom panel        ``'bottom'``, ``'b'``
            top panel           ``'top'``, ``'t'``
            "best" inset        ``'best'``, ``'inset'``, ``'i'``, ``0``
            upper right inset   ``'upper right'``, ``'ur'``, ``1``
            upper left inset    ``'upper left'``, ``'ul'``, ``2``
            lower left inset    ``'lower left'``, ``'ll'``, ``3``
            lower right inset   ``'lower right'``, ``'lr'``, ``4``
            center left inset   ``'center left'``, ``'cl'``, ``5``
            center right inset  ``'center right'``, ``'cr'``, ``6``
            lower center inset  ``'lower center'``, ``'lc'``, ``7``
            upper center inset  ``'upper center'``, ``'uc'``, ``8``
            center inset        ``'center'``, ``'c'``, ``9``
            ==================  =======================================

        width : float or str, optional
            The space allocated for outer legends. This does nothing
            if :rcraw:`tight` is ``True``. Units are interpreted by
            `~proplot.utils.units`.
        space : float or str, optional
            The space between the axes and the legend for outer legends.
            Units are interpreted by `~proplot.utils.units`.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, defaut is :rc:`subplots.panelspace`.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~proplot.wrappers.legend_wrapper`.
        """
        loc = self._loc_translate(loc, width=width, space=space)
        if isinstance(loc, np.ndarray):
            loc = loc.tolist()

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
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
            frameon = _notNone(kwargs.get('frame', None), kwargs.get(
                'frameon', None), rc['legend.frameon'])
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
        return legend_wrapper(self, *args, loc=loc, **kwargs)

    def draw(self, renderer=None, *args, **kwargs):
        """Adds post-processing steps before axes is drawn."""
        self._reassign_title()
        super().draw(renderer, *args, **kwargs)

    def get_size_inches(self):
        """Returns the width and the height of the axes in inches."""
        width, height = self.figure.get_size_inches()
        width = width * abs(self.get_position().width)
        height = height * abs(self.get_position().height)
        return width, height

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Adds post-processing steps before tight bounding box is
        calculated, and stores the bounding box as an attribute."""
        self._reassign_title()
        bbox = super().get_tightbbox(renderer, *args, **kwargs)
        self._tightbbox = bbox
        return bbox

    def heatmap(self, *args, **kwargs):
        """Calls `~matplotlib.axes.Axes.pcolormesh` and applies default formatting
        that is suitable for heatmaps: no gridlines, no minor ticks, and major
        ticks at the center of each grid box."""
        obj = self.pcolormesh(*args, **kwargs)
        xlocator, ylocator = None, None
        if hasattr(obj, '_coordinates'):
            coords = obj._coordinates
            coords = (coords[1:, ...] + coords[:-1, ...]) / 2
            coords = (coords[:, 1:, :] + coords[:, :-1, :]) / 2
            xlocator, ylocator = coords[0, :, 0], coords[:, 0, 1]
        self.format(
            xgrid=False, ygrid=False, xtickminor=False, ytickminor=False,
            xlocator=xlocator, ylocator=ylocator,
        )
        return obj

    def inset_axes(self, bounds, *, transform=None, zorder=4,
                   zoom=True, zoom_kw=None, **kwargs):
        """
        Like the builtin `~matplotlib.axes.Axes.inset_axes` method, but
        draws an inset `XYAxes` axes and adds some options.

        Parameters
        ----------
        bounds : list of float
            The bounds for the inset axes, listed as ``(x, y, width, height)``.
        transform : {'data', 'axes', 'figure'} or \
`~matplotlib.transforms.Transform`, optional
            The transform used to interpret `bounds`. Can be a
            `~matplotlib.transforms.Transform` object or a string representing
            the `~matplotlib.axes.Axes.transData`,
            `~matplotlib.axes.Axes.transAxes`,
            or `~matplotlib.figure.Figure.transFigure` transforms. Default is
            ``'axes'``, i.e. `bounds` is in axes-relative coordinates.
        zorder : float, optional
            The `zorder \
<https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
            of the axes, should be greater than the zorder of
            elements in the parent axes. Default is ``4``.
        zoom : bool, optional
            Whether to draw lines indicating the inset zoom using
            `~Axes.indicate_inset_zoom`. The lines will automatically
            adjust whenever the parent axes or inset axes limits are changed.
            Default is ``True``.
        zoom_kw : dict, optional
            Passed to `~Axes.indicate_inset_zoom`.

        Other parameters
        ----------------
        **kwargs
            Passed to `XYAxes`.
        """
        # Carbon copy with my custom axes
        if not transform:
            transform = self.transAxes
        else:
            transform = _get_transform(self, transform)
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        locator = self._make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = XYAxes(self.figure, bb.bounds,
                    zorder=zorder, label=label, **kwargs)
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
        linewidth = _notNone(
            lw, linewidth, rc['axes.linewidth'],
            names=('lw', 'linewidth'))
        edgecolor = _notNone(
            color, edgecolor, rc['axes.edgecolor'],
            names=('color', 'edgecolor'))
        if not parent:
            raise ValueError(f'{self} is not an inset axes.')
        xlim, ylim = self.get_xlim(), self.get_ylim()
        rect = (xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0])

        # Call indicate_inset
        rectpatch, connects = parent.indicate_inset(
            rect, self, linewidth=linewidth, edgecolor=edgecolor, alpha=alpha,
            **kwargs)

        # Update zoom or adopt properties from old one
        if self._inset_zoom_data:
            rectpatch_old, connects_old = self._inset_zoom_data
            rectpatch.update_from(rectpatch_old)
            rectpatch_old.set_visible(False)
            for line, line_old in zip(connects, connects_old):
                visible = line.get_visible()
                line.update_from(line_old)
                line.set_linewidth(line_old.get_linewidth())
                line.set_visible(visible)
                line_old.set_visible(False)
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
            The panel width. Units are interpreted by `~proplot.utils.units`.
            Default is :rc:`subplots.panelwidth`.
        space : float or str or list thereof, optional
            Empty space between the main subplot and the panel.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, defaut is :rc:`subplots.panelspace`.
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

    @_standardize_1d
    @_cmap_changer
    def parametric(self, *args, values=None,
                   cmap=None, norm=None,
                   interp=0, **kwargs):
        """
        Invoked when you pass the `cmap` keyword argument to
        `~matplotlib.axes.Axes.plot`. Draws a "colormap line",
        i.e. a line whose color changes as a function of the parametric
        coordinate ``values``. using the input colormap ``cmap``.

        This is actually a collection of lines, added as a
        `~matplotlib.collections.LineCollection` instance. See
        `this matplotlib example \
<https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line>`__.

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
            If greater than ``0``, we interpolate to additional points
            between the `values` coordinates. The number corresponds to the
            number of additional color levels between the line joints
            and the halfway points between line joints.
        """
        # First error check
        # WARNING: So far this only works for 1D *x* and *y* coordinates.
        # Cannot draw multiple colormap lines at once
        if values is None:
            raise ValueError('Requires a "values" keyword arg.')
        if len(args) not in (1, 2):
            raise ValueError(f'Requires 1-2 arguments, got {len(args)}.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(
            y.shape[-1]) if len(args) == 1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim != 1 or y.ndim != 1 or values.ndim != 1:
            raise ValueError(
                f'x ({x.ndim}d), y ({y.ndim}d), and values ({values.ndim}d)'
                ' must be 1-dimensional.')
        if len(x) != len(y) or len(x) != len(values) or len(y) != len(values):
            raise ValueError(
                f'{len(x)} xs, {len(y)} ys, but {len(values)} '
                ' colormap values.')

        # Interpolate values to allow for smooth gradations between values
        # (bins=False) or color switchover halfway between points (bins=True)
        # Then optionally interpolate the corresponding colormap values
        if interp > 0:
            xorig, yorig, vorig = x, y, values
            x, y, values = [], [], []
            for j in range(xorig.shape[0] - 1):
                idx = (
                    slice(None, -1) if j + 1 < xorig.shape[0] - 1
                    else slice(None))
                x.extend(np.linspace(
                    xorig[j], xorig[j + 1], interp + 2)[idx].flat)
                y.extend(np.linspace(
                    yorig[j], yorig[j + 1], interp + 2)[idx].flat)
                values.extend(np.linspace(
                    vorig[j], vorig[j + 1], interp + 2)[idx].flat)
            x, y, values = np.array(x), np.array(y), np.array(values)
        coords = []
        levels = utils.edges(values)
        for j in range(y.shape[0]):
            # Get x/y coordinates and values for points to the 'left' and
            # 'right' of each joint
            if j == 0:
                xleft, yleft = [], []
            else:
                xleft = [(x[j - 1] + x[j]) / 2, x[j]]
                yleft = [(y[j - 1] + y[j]) / 2, y[j]]
            if j + 1 == y.shape[0]:
                xright, yright = [], []
            else:
                xleft = xleft[:-1]  # prevent repetition when joined with right
                yleft = yleft[:-1]
                xright = [x[j], (x[j + 1] + x[j]) / 2]
                yright = [y[j], (y[j + 1] + y[j]) / 2]
            pleft = np.stack((xleft, yleft), axis=1)
            pright = np.stack((xright, yright), axis=1)
            coords.append(np.concatenate((pleft, pright), axis=0))

        # Create LineCollection and update with values
        hs = mcollections.LineCollection(
            np.array(coords), cmap=cmap, norm=norm,
            linestyles='-', capstyle='butt', joinstyle='miter')
        hs.set_array(np.array(values))
        hs.update({key: value for key, value in kwargs.items()
                   if key not in ('color',)})

        # Add collection, with some custom attributes
        self.add_collection(hs)
        if self.get_autoscale_on() and self.ignore_existing_data_limits:
            self.autoscale_view()  # data limits not updated otherwise
        hs.values = values
        hs.levels = levels  # needed for other functions some
        return hs

    def violins(self, *args, **kwargs):
        """Alias for `~matplotlib.axes.Axes.violinplot`."""
        return self.violinplot(*args, **kwargs)

    #: Alias for `~Axes.panel_axes`.
    panel = panel_axes

    #: Alias for `~Axes.inset_axes`.
    inset = inset_axes

    @property
    def number(self):
        """The axes number, controls a-b-c label order and order of
        appearence in the `~proplot.subplots.subplot_grid` returned by
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

    # Wrapped by special functions
    # Also support redirecting to Basemap methods
    text = _text_wrapper(
        maxes.Axes.text
    )
    plot = _plot_wrapper(_standardize_1d(_add_errorbars(_cycle_changer(
        maxes.Axes.plot
    ))))
    scatter = _scatter_wrapper(_standardize_1d(_add_errorbars(_cycle_changer(
        maxes.Axes.scatter
    ))))
    bar = _bar_wrapper(_standardize_1d(_add_errorbars(_cycle_changer(
        maxes.Axes.bar
    ))))
    barh = _barh_wrapper(
        maxes.Axes.barh
    )  # calls self.bar
    hist = _hist_wrapper(_standardize_1d(_cycle_changer(
        maxes.Axes.hist
    )))
    boxplot = _boxplot_wrapper(_standardize_1d(_cycle_changer(
        maxes.Axes.boxplot
    )))
    violinplot = _violinplot_wrapper(_standardize_1d(_add_errorbars(
        _cycle_changer(maxes.Axes.violinplot)
    )))
    fill_between = _fill_between_wrapper(_standardize_1d(_cycle_changer(
        maxes.Axes.fill_between
    )))
    fill_betweenx = _fill_betweenx_wrapper(_standardize_1d(_cycle_changer(
        maxes.Axes.fill_betweenx
    )))

    # Wrapped by cycle wrapper and standardized
    pie = _standardize_1d(_cycle_changer(
        maxes.Axes.pie
    ))
    stem = _standardize_1d(_cycle_changer(
        maxes.Axes.stem
    ))
    step = _standardize_1d(_cycle_changer(
        maxes.Axes.step
    ))

    # Wrapped by cmap wrapper and standardized
    # Also support redirecting to Basemap methods
    hexbin = _standardize_1d(_cmap_changer(
        maxes.Axes.hexbin
    ))
    contour = _standardize_2d(_cmap_changer(
        maxes.Axes.contour
    ))
    contourf = _standardize_2d(_cmap_changer(
        maxes.Axes.contourf
    ))
    pcolor = _standardize_2d(_cmap_changer(
        maxes.Axes.pcolor
    ))
    pcolormesh = _standardize_2d(_cmap_changer(
        maxes.Axes.pcolormesh
    ))
    quiver = _standardize_2d(_cmap_changer(
        maxes.Axes.quiver
    ))
    streamplot = _standardize_2d(_cmap_changer(
        maxes.Axes.streamplot
    ))
    barbs = _standardize_2d(_cmap_changer(
        maxes.Axes.barbs
    ))
    imshow = _cmap_changer(
        maxes.Axes.imshow
    )

    # Wrapped only by cmap wrapper
    tripcolor = _cmap_changer(
        maxes.Axes.tripcolor
    )
    tricontour = _cmap_changer(
        maxes.Axes.tricontour
    )
    tricontourf = _cmap_changer(
        maxes.Axes.tricontourf
    )
    hist2d = _cmap_changer(
        maxes.Axes.hist2d
    )
    spy = _cmap_changer(
        maxes.Axes.spy
    )
    matshow = _cmap_changer(
        maxes.Axes.matshow
    )


dualxy_kwargs = (
    'label', 'locator', 'formatter', 'ticks', 'ticklabels',
    'minorlocator', 'minorticks', 'tickminor',
    'ticklen', 'tickrange', 'tickdir', 'ticklabeldir', 'tickrotation',
    'bounds', 'margin', 'color', 'grid', 'gridminor', 'gridcolor',
)

dualxy_descrip = """
Makes a secondary *%(x)s* axis for denoting equivalent *%(x)s*
coordinates in *alternate units*.

Parameters
----------
forward : function, optional
    Function used to transform units from the original axis to the
    secondary axis. Should take 1 value and perform a *forward
    linear transformation*. For example, to convert Kelvin to Celsius,
    use ``ax.dual%(x)s(lambda x: x - 273.15)``. To convert kilometers to
    meters, use ``ax.dual%(x)s(lambda x: x*1e3)``.
inverse : function, optional
    Function used to transform units from the secondary axis back to
    the original axis. If `forward` was a non-linear function, you
    *must* provide this, or the transformation will be incorrect!

    For example, to apply the square, use
    ``ax.dual%(x)s(lambda x: x**2, lambda x: x**0.5)``.
scale : scale-spec, optional
    The axis scale from which forward and inverse transformations
    are inferred. Passed to `~proplot.axistools.Scale`.

    For example, to apply the inverse, use ``ax.dual%(x)s('inverse')``;
    To apply the base-10 exponential function, use
    ``ax.dual%(x)s(('exp', 10, 1, 10))``
    or ``ax.dual%(x)s(plot.Scale('exp', 10))``.
scale_kw : dict-like, optional
    Ignored if `scale` is ``None``. Passed to
    `~proplot.axistools.Scale`.
%(args)s : optional
    Prepended with ``'y'`` and passed to `Axes.format`.
"""

altxy_descrip = """
Alias and more intuitive name for `~XYAxes.twin%(y)s`.
The matplotlib `~matplotlib.axes.Axes.twiny` function
generates two *x* axes with a shared ("twin") *y* axis.
Enforces the following settings.

* Places the old *%(x)s* axis on the %(x1)s and the new *%(x)s* axis
  on the %(x2)s.
* Makes the old %(x2)s spine invisible and the new %(x1)s, %(y1)s,
  and %(y2)s spines invisible.
* Adjusts the *%(x)s* axis tick, tick label, and axis label positions
  according to the visible spine positions.
* Locks the old and new *%(y)s* axis limits and scales, and makes the new
  %(y)s axis labels invisible.
"""

twinxy_descrip = """
Mimics matplotlib's `~matplotlib.axes.Axes.twin%(y)s`.
Enforces the following settings.

* Places the old *%(x)s* axis on the %(x1)s and the new *%(x)s* axis
  on the %(x2)s.
* Makes the old %(x2)s spine invisible and the new %(x1)s, %(y1)s,
  and %(y2)s spines invisible.
* Adjusts the *%(x)s* axis tick, tick label, and axis label positions
  according to the visible spine positions.
* Locks the old and new *%(y)s* axis limits and scales, and makes the new
  %(y)s axis labels invisible.
"""


def _parse_dualxy_args(x, transform, transform_kw, kwargs):
    """Interprets the dualx and dualy transform and various keyword
    arguments. Returns a list of forward transform, inverse transform, and
    overrides for default locators and formatters."""
    # Transform using input functions
    # TODO: Also support transforms? Probably not -- transforms are a huge
    # group that include ND and non-invertable transformations, but transforms
    # used for axis scales are subset of invertible 1D functions
    funcscale_kw = {}
    transform_kw = transform_kw or {}
    if isinstance(transform, (str, mscale.ScaleBase)) or transform_kw:
        transform = transform or 'linear'
        scale = axistools.Scale(transform, **transform_kw)
        transform = scale.get_transform()
        funcscale_funcs = (transform.transform, transform.inverted().transform)
        for key in ('major_locator', 'minor_locator',
                    'major_formatter', 'minor_formatter'):
            default = getattr(scale, '_' + key, None)
            if default:
                funcscale_kw[key] = default
    elif (np.iterable(transform) and len(transform) == 2
          and all(callable(itransform) for itransform in transform)):
        funcscale_funcs = transform
    elif callable(transform):
        funcscale_funcs = (transform, lambda x: x)
    else:
        raise ValueError(
            f'Invalid transform {transform!r}. '
            'Must be function, tuple of two functions, or scale name.')
    # Parse keyword args intended for format() command
    kwargs_bad = {}
    for key in (*kwargs.keys(),):
        value = kwargs.pop(key)
        if key[0] == x and key[1:] in dualxy_kwargs:
            warnings.warn(
                f'dual{x}() keyword arg {key!r} is deprecated. '
                f'Use {key[1:]!r} instead.')
            kwargs[key] = value
        elif key in dualxy_kwargs:
            kwargs[x + key] = value
        else:
            kwargs_bad[key] = value
        if kwargs_bad:
            raise TypeError(
                f'dual{x}() got unexpected keyword argument(s): {kwargs_bad}')
    return funcscale_funcs, funcscale_kw, kwargs


def _rcloc_to_stringloc(x, string):  # figures out string location
    """Gets *location string* from the *boolean* "left", "right", "top", and
    "bottom" rc settings, e.g. :rc:`axes.spines.left` or :rc:`ytick.left`.
    Might be ``None`` if settings are unchanged."""
    # For x axes
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
    # For y axes
    else:
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


class XYAxes(Axes):
    """
    Axes subclass for ordinary 2D cartesian coordinates. Adds several new
    methods and overrides existing ones.
    """
    #: The registered projection name.
    name = 'xy'

    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`, `Axes`
        """
        # Impose default formatter
        super().__init__(*args, **kwargs)
        formatter = axistools.Formatter('auto')
        self.xaxis.set_major_formatter(formatter)
        self.yaxis.set_major_formatter(formatter)
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True
        # Custom attributes
        self._datex_rotated = False  # whether manual rotation has been applied
        self._dualy_data = None  # for scaling units on opposite side of ax
        self._dualx_data = None

    def _altx_overrides(self):
        """Applies alternate *x* axis overrides."""
        # Unlike matplotlib API, we strong arm user into certain twin axes
        # settings... doesn't really make sense to have twin axes without this
        if self._altx_child is not None:  # altx was called on this axes
            self._shared_y_axes.join(self, self._altx_child)
            self.spines['top'].set_visible(False)
            self.spines['bottom'].set_visible(True)
            self.xaxis.tick_bottom()
            self.xaxis.set_label_position('bottom')
        elif self._altx_parent is not None:  # this axes is the result of altx
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
            self._shared_x_axes.join(self, self._alty_child)
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
        if (not isinstance(self.xaxis.converter, mdates.DateConverter)
                or self._datex_rotated):
            return
        rotation = rc['axes.formatter.timerotation']
        kw = {'rotation': rotation}
        if rotation not in (0, 90, -90):
            kw['ha'] = ('right' if rotation > 0 else 'left')
        for label in self.xaxis.get_ticklabels():
            label.update(kw)
        self._datex_rotated = True  # do not need to apply more than once

    def _dualx_overrides(self):
        """Locks child "dual" *x* axis limits to the parent."""
        # Why did I copy and paste the dualx/dualy code you ask? Copy
        # pasting is bad, but so are a bunch of ugly getattr(attr)() calls
        data = self._dualx_data
        if data is None:
            return
        funcscale_funcs, funcscale_kw = data
        # Build the FuncScale
        # Sometimes we do *not* want to apply default locator and formatter
        # overrides, in case user has manually changed them! Also, sometimes
        # we want to borrow method that sets default from the scale from which
        # forward and inverse funcs were drawn, instead of from FuncScale.
        child = self._altx_child
        scale_parent = self.xaxis._scale
        scale_func = axistools.Scale(
            'function', funcscale_funcs[::-1],
            scale_parent.get_transform(),
            **funcscale_kw
        )
        scale_defaults = scale_func if isinstance(
            scale_parent, mscale.LinearScale) else scale_parent
        try:
            scale_defaults.set_default_locators_and_formatters(
                child.xaxis, only_if_default=True,
            )
        except TypeError:
            pass  # do nothing if axis has native matplotlib scale
        child.xaxis._scale = scale_func
        child._update_transScale()
        child.stale = True
        child.autoscale_view(scaley=False)
        # Transform axis limits
        # If the transform flipped the limits, when we set axis limits, it
        # will get flipped again! So reverse the flip
        lim = self.get_xlim()
        nlim = list(map(funcscale_funcs[0], np.array(lim)))
        if np.sign(np.diff(lim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]
        child.set_xlim(nlim)

    def _dualy_overrides(self):
        """Locks child "dual" *y* axis limits to the parent."""
        data = self._dualy_data
        if data is None:
            return
        funcscale_funcs, funcscale_kw = data
        child = self._alty_child
        scale_parent = self.yaxis._scale
        scale_func = axistools.Scale(
            'function', funcscale_funcs[::-1], scale_parent.get_transform(),
            **funcscale_kw
        )
        scale_defaults = scale_func if isinstance(
            scale_parent, mscale.LinearScale) else scale_parent
        try:
            scale_defaults.set_default_locators_and_formatters(
                child.xaxis, only_if_default=True,
            )
        except TypeError:
            pass
        child.yaxis._scale = scale_func
        child._update_transScale()
        child.stale = True
        child.autoscale_view(scalex=False)
        lim = self.get_ylim()
        nlim = list(map(funcscale_funcs[0], np.array(lim)))
        if np.sign(np.diff(lim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]
        child.set_ylim(nlim)

    def _hide_labels(self):
        """Function called at drawtime that enforces "shared" axis and
        tick labels. If this is not called at drawtime, "shared" labels can
        be inadvertantly turned off e.g. when the axis scale is changed."""
        for x in 'xy':
            # "Shared" axis and tick labels
            axis = getattr(self, x + 'axis')
            share = getattr(self, '_share' + x)
            if share is not None:
                level = getattr(self, '_share' + x + '_level')
                if level > 0:
                    axis.label.set_visible(False)
                if level > 2:
                    axis.set_major_formatter(mticker.NullFormatter())
            # Enforce no minor ticks labels
            # TODO: Document this?
            axis.set_minor_formatter(mticker.NullFormatter())

    def _sharex_setup(self, sharex, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        # Call Axes method
        super()._sharex_setup(sharex, level)  # sets up panels
        if sharex in (None, self) or not isinstance(sharex, XYAxes):
            return
        # Builtin sharing features
        if level > 0:
            self._sharex = sharex
        if level > 1:
            self._shared_x_axes.join(self, sharex)

    def _sharey_setup(self, sharey, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        # Call Axes method
        super()._sharey_setup(sharey, level)
        if sharey in (None, self) or not isinstance(sharey, XYAxes):
            return
        # Builtin features
        if level > 0:
            self._sharey = sharey
        if level > 1:
            self._shared_y_axes.join(self, sharey)

    def format(
            self, *,
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
            xreverse=None, yreverse=None,
            xlabel=None, ylabel=None,
            xlim=None, ylim=None,
            xscale=None, yscale=None,
            xrotation=None, yrotation=None,
            xformatter=None, yformatter=None,
            xticklabels=None, yticklabels=None,
            xticks=None, xminorticks=None,
            xlocator=None, xminorlocator=None,
            yticks=None, yminorticks=None,
            ylocator=None, yminorlocator=None,
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
            `~matplotlib.axes.Axes.set_xlim` and
            `~matplotlib.axes.Axes.set_ylim`.
        xreverse, yreverse : bool, optional
            Sets whether the *x* and *y* axis are oriented in the "reverse"
            direction. The "normal" direction is increasing to the right for
            the *x* axis and to the top for the *y* axis. The "reverse"
            direction is increasing to the left for the *x* axis and to the
            bottom for the *y* axis.
        xscale, yscale : axis scale spec, optional
            The *x* and *y* axis scales. Passed to the
            `~proplot.axistools.Scale` constructor. For example,
            ``xscale='log'`` applies logarithmic scaling, and
            ``xscale=('cutoff', 0.5, 2)`` applies scaling according to the
            class generated by `plot.CutoffScaleFactory('cutoff', 0.5, 2)`.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.axistools.Scale`.
        xspineloc, yspineloc : {'both', 'bottom', 'top', 'left', 'right', \
'neither', 'center', 'zero'}, optional
            The *x* and *y* axis spine locations.
        xloc, yloc : optional
            Aliases for `xspineloc`, `yspineloc`.
        xtickloc, ytickloc : {'both', 'bottom', 'top', 'left', 'right', \
'neither'}, optional
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
            The rotation for *x* and *y* axis tick labels. Default is ``0``
            for normal axes, :rc:`axes.formatter.timerotation` for time
            *x* axes.
        xtickrange, ytickrange : (float, float), optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``(-1,1)`` with
            axis range ``(-5,5)`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1.
        xmargin, ymargin : float, optional
            The default margin between plotted content and the *x* and *y* axis
            spines. Value is proportional to the width, height of the axes.
            Use this if you want whitespace between plotted content
            and the spines, but don't want to explicitly set `xlim` or `ylim`.
        xbounds, ybounds : (float, float), optional
            The *x* and *y* axis data bounds within which to draw the spines.
            For example, the axis range ``(0, 4)`` with bounds ``(1, 4)``
            will prevent the spines from meeting at the origin.
        xcolor, ycolor : color-spec, optional
            Color for the *x* and *y* axis spines, ticks, tick labels, and axis
            labels. Default is :rc:`color`. Use e.g. ``ax.format(color='red')``
            to set for both axes.
        xticklen, yticklen : float or str, optional
            Tick lengths for the *x* and *y* axis. Units are interpreted by
            `~proplot.utils.units`, with "points" as the numeric unit. Default
            is :rc:`ticklen`.

            Minor tick lengths are scaled according
            to :rc:`ticklenratio`. Use e.g. ``ax.format(ticklen=1)`` to
            set for both axes.
        fixticks : bool, optional
            Whether to always transform the tick locators to a
            `~matplotlib.ticker.FixedLocator` instance. Default is ``False``.
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
        If you plot something with a `datetime64 \
<https://docs.scipy.org/doc/numpy/reference/arrays.datetime.html>`__,
        `pandas.Timestamp`, `pandas.DatetimeIndex`, `datetime.date`,
        `datetime.time`, or `datetime.datetime` array as the *x* or *y* axis
        coordinate, the axis ticks and tick labels will be automatically
        formatted as dates.

        See also
        --------
        :py:obj:`Axes.format`, :py:obj:`Axes.context`
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
            xlabel_kw = xlabel_kw or {}
            ylabel_kw = ylabel_kw or {}
            xscale_kw = xscale_kw or {}
            yscale_kw = yscale_kw or {}
            xlocator_kw = xlocator_kw or {}
            ylocator_kw = ylocator_kw or {}
            xformatter_kw = xformatter_kw or {}
            yformatter_kw = yformatter_kw or {}
            xminorlocator_kw = xminorlocator_kw or {}
            yminorlocator_kw = yminorlocator_kw or {}
            # Flexible keyword args, declare defaults
            xmargin = _notNone(xmargin, rc['axes.xmargin'])
            ymargin = _notNone(ymargin, rc['axes.ymargin'])
            xtickdir = _notNone(xtickdir, rc['xtick.direction'])
            ytickdir = _notNone(ytickdir, rc['ytick.direction'])
            xtickminor = _notNone(xtickminor, rc['xtick.minor.visible'])
            ytickminor = _notNone(ytickminor, rc['ytick.minor.visible'])
            xformatter = _notNone(xticklabels, xformatter,
                                  None, names=('xticklabels', 'xformatter'))
            yformatter = _notNone(yticklabels, yformatter,
                                  None, names=('yticklabels', 'yformatter'))
            xlocator = _notNone(xticks, xlocator, None,
                                names=('xticks', 'xlocator'))
            ylocator = _notNone(yticks, ylocator, None,
                                names=('yticks', 'ylocator'))
            xminorlocator = _notNone(xminorticks, xminorlocator, None, names=(
                'xminorticks', 'xminorlocator'))
            yminorlocator = _notNone(yminorticks, yminorlocator, None, names=(
                'yminorticks', 'yminorlocator'))
            # Grid defaults are more complicated
            axis = rc.get('axes.grid.axis')  # always need this property
            grid, which = rc['axes.grid'], rc['axes.grid.which']
            if which is not None or grid is not None:  # if *one* was changed
                if grid is None:
                    grid = rc.get('axes.grid')
                elif which is None:
                    which = rc.get('axes.grid.which')
                xgrid = _notNone(xgrid,
                                 grid and axis in ('x', 'both')
                                 and which in ('major', 'both'))
                ygrid = _notNone(ygrid,
                                 grid and axis in ('y', 'both')
                                 and which in ('major', 'both'))
                xgridminor = _notNone(xgridminor, grid
                                      and axis in ('x', 'both')
                                      and which in ('minor', 'both'))
                ygridminor = _notNone(ygridminor, grid
                                      and axis in ('y', 'both')
                                      and which in ('minor', 'both'))

            # Sensible defaults for spine, tick, tick label, and label locs
            # NOTE: Allow tick labels to be present without ticks! User may
            # want this sometimes! Same goes for spines!
            xspineloc = _notNone(
                xloc, xspineloc, None, names=('xloc', 'xspineloc')
            )
            yspineloc = _notNone(
                yloc, yspineloc, None, names=('yloc', 'yspineloc')
            )
            xtickloc = _notNone(
                xtickloc, xspineloc, _rcloc_to_stringloc('x', 'xtick')
            )
            ytickloc = _notNone(
                ytickloc, yspineloc, _rcloc_to_stringloc('y', 'ytick')
            )
            xspineloc = _notNone(
                xspineloc, _rcloc_to_stringloc('x', 'axes.spines')
            )
            yspineloc = _notNone(
                yspineloc, _rcloc_to_stringloc('y', 'axes.spines')
            )
            if xtickloc != 'both':
                xticklabelloc = _notNone(xticklabelloc, xtickloc)
                xlabelloc = _notNone(xlabelloc, xticklabelloc)
                if xlabelloc not in (None, 'bottom', 'top'):  # e.g. "both"
                    xlabelloc = 'bottom'
            if ytickloc != 'both':
                yticklabelloc = _notNone(yticklabelloc, ytickloc)
                ylabelloc = _notNone(ylabelloc, yticklabelloc)
                if ylabelloc not in (None, 'left', 'right'):
                    ylabelloc = 'left'

            # Begin loop
            for (
                x, axis,
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
                ('x', 'y'), (self.xaxis, self.yaxis),
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
                (xlocator_kw, ylocator_kw),
                (xminorlocator_kw, yminorlocator_kw),
                (xformatter_kw, yformatter_kw),
            ):
                # Axis limits
                # NOTE: 3.1+ has axis.set_inverted(), below is from source code
                if lim is not None:
                    getattr(self, 'set_' + x + 'lim')(lim)
                if reverse is not None:
                    lo, hi = axis.get_view_interval()
                    if reverse:
                        lim = (max(lo, hi), min(lo, hi))
                    else:
                        lim = (min(lo, hi), max(lo, hi))
                    axis.set_view_interval(*lim, ignore=True)
                # Axis scale
                # WARNING: This relies on monkey patch of mscale.scale_factory
                # that allows it to accept a custom scale class!
                # WARNING: Changing axis scale also changes default locators
                # and formatters, so do it first
                if scale is not None:
                    scale = axistools.Scale(scale, **scale_kw)
                    getattr(self, 'set_' + x + 'scale')(scale)
                # Is this a date axis?
                # NOTE: Make sure to get this *after* lims set!
                # See: https://matplotlib.org/api/units_api.html
                # And: https://matplotlib.org/api/dates_api.html
                # Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axis.py # noqa
                # The axis_date() method just applies DateConverter
                date = isinstance(axis.converter, mdates.DateConverter)

                # Fix spines
                kw = rc.fill({
                    'linewidth': 'axes.linewidth',
                    'color': 'axes.edgecolor',
                })
                if color is not None:
                    kw['color'] = color
                sides = ('bottom', 'top') if x == 'x' else ('left', 'right')
                spines = [self.spines[s] for s in sides]
                for spine, side in zip(spines, sides):
                    # Line properties
                    # Override if we're settings spine bounds
                    # In this case just have spines on edges by default
                    if bounds is not None and spineloc not in sides:
                        spineloc = sides[0]
                    # Eliminate sides
                    if spineloc == 'neither':
                        spine.set_visible(False)
                    elif spineloc == 'both':
                        spine.set_visible(True)
                    elif spineloc in sides:  # make relevant spine visible
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
                                raise ValueError(
                                    f'Invalid {x} spine location {spineloc!r}.'
                                    f' Options are '
                                    + ', '.join(map(
                                        repr, (*sides, 'both', 'neither')
                                    )) + '.'
                                )
                    # Apply spine bounds
                    if bounds is not None and spine.get_visible():
                        spine.set_bounds(*bounds)
                    spine.update(kw)
                # Get available spines, needed for setting tick locations
                spines = [side for side, spine in zip(
                    sides, spines) if spine.get_visible()]

                # Tick and grid settings for major and minor ticks separately
                # Override is just a "new default", but user can override this
                def grid_dict(grid):
                    return {
                        'grid_color': grid + '.color',
                        'grid_alpha': grid + '.alpha',
                        'grid_linewidth': grid + '.linewidth',
                        'grid_linestyle': grid + '.linestyle',
                    }
                for which, igrid in zip(('major', 'minor'), (grid, gridminor)):
                    # Tick properties
                    kw_ticks = rc.category(x + 'tick.' + which)
                    if kw_ticks is None:
                        kw_ticks = {}
                    else:
                        kw_ticks.pop('visible', None)  # invalid setting
                    if ticklen is not None:
                        if which == 'major':
                            kw_ticks['size'] = utils.units(ticklen, 'pt')
                        else:
                            kw_ticks['size'] = utils.units(
                                ticklen, 'pt') * rc.get('ticklenratio')
                    # Grid style and toggling
                    if igrid is not None:
                        # toggle with special global props
                        axis.grid(igrid, which=which)
                    if which == 'major':
                        kw_grid = rc.fill(grid_dict('grid'))
                    else:
                        kw_major = kw_grid
                        kw_grid = rc.fill(grid_dict('gridminor'))
                        kw_grid.update({
                            key: value for key, value in kw_major.items()
                            if key not in kw_grid
                        })
                    # Changed rc settings
                    axis.set_tick_params(which=which, **kw_grid, **kw_ticks)

                # Tick and ticklabel properties that apply to major and minor
                # * Weird issue causes set_tick_params to reset/forget grid
                #   is turned on if you access tick.gridOn directly, instead of
                #   passing through tick_params. Since gridOn is undocumented
                #   feature, don't use it. So calling _format_axes() a second
                #   time will remove the lines.
                # * Can specify whether the left/right/bottom/top spines get
                #   ticks; sides will be group of left/right or top/bottom.
                # * Includes option to draw spines but not draw ticks on that
                #   spine, e.g. on the left/right edges
                kw = {}
                translate = {None: None, 'both': sides,
                             'neither': (), 'none': ()}
                if bounds is not None and tickloc not in sides:
                    tickloc = sides[0]  # override to just one side
                ticklocs = translate.get(tickloc, (tickloc,))
                if ticklocs is not None:
                    kw.update({side: (side in ticklocs) for side in sides})
                kw.update({  # override
                    side: False for side in sides if side not in spines
                })
                # Tick label sides
                # Will override to make sure only appear where ticks are
                ticklabellocs = translate.get(ticklabelloc, (ticklabelloc,))
                if ticklabellocs is not None:
                    kw.update({
                        'label' + side: (side in ticklabellocs)
                        for side in sides
                    })
                kw.update({  # override
                    'label' + side: False for side in sides
                    if (side not in spines or (ticklocs is not None
                                               and side not in ticklocs))
                })  # override
                # The axis label side
                if labelloc is None:
                    if ticklocs is not None:
                        options = [side for side in sides if (
                            side in ticklocs and side in spines)]
                        if len(options) == 1:
                            labelloc = options[0]
                elif labelloc not in sides:
                    raise ValueError(
                        f'Got labelloc {labelloc!r}, valid options are '
                        + ', '.join(map(repr, sides)) + '.'
                    )
                # Apply
                axis.set_tick_params(which='both', **kw)
                if labelloc is not None:
                    axis.set_label_position(labelloc)

                # Tick label settings
                # First color and size
                kw = rc.fill({
                    'labelcolor': 'tick.labelcolor',  # new props
                    'labelsize': 'tick.labelsize',
                    'color': x + 'tick.color',
                })
                if color:
                    kw['color'] = color
                    kw['labelcolor'] = color
                # Tick direction and rotation
                if tickdir == 'in':
                    kw['pad'] = 1  # ticklabels should be much closer
                if ticklabeldir == 'in':  # put tick labels inside the plot
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
                    'weight': 'tick.labelweight'
                })
                if rotation is not None:
                    kw = {'rotation': rotation}
                    if x == 'x':
                        self._datex_rotated = True
                        if rotation not in (0, 90, -90):
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
                    'color': 'axes.edgecolor',
                    'weight': 'axes.labelweight',
                    'fontsize': 'axes.labelsize',
                    'fontfamily': 'font.family',
                })
                if label is not None:
                    kw['text'] = label
                if color:
                    kw['color'] = color
                kw.update(label_kw)
                if kw:  # NOTE: initially keep spanning labels off
                    self._update_axislabels(x, **kw)

                # Major and minor locator
                # WARNING: MultipleLocator fails sometimes, notably when doing
                # boxplot. Tick labels moved to left and are incorrect.
                if locator is not None:
                    locator = axistools.Locator(locator, **locator_kw)
                    axis.set_major_locator(locator)
                    if isinstance(locator, mticker.IndexLocator):
                        tickminor = False  # 'index' minor ticks make no sense
                if not tickminor and tickminorlocator is None:
                    axis.set_minor_locator(axistools.Locator('null'))
                elif tickminorlocator is not None:
                    axis.set_minor_locator(axistools.Locator(
                        tickminorlocator, **minorlocator_kw))

                # Major formatter
                # NOTE: Only reliable way to disable ticks labels and then
                # restore them is by messing with the formatter, *not* setting
                # labelleft=False, labelright=False, etc.
                if (formatter is not None or tickrange is not None) and not (
                        isinstance(axis.get_major_formatter(),
                                   mticker.NullFormatter)
                        and getattr(self, '_share' + x)):
                    # Tick range
                    if tickrange is not None:
                        if formatter not in (None, 'auto'):
                            warnings.warn(
                                'The tickrange feature requires '
                                'proplot.AutoFormatter formatter. Overriding '
                                'input formatter.'
                            )
                        formatter = 'auto'
                        formatter_kw.setdefault('tickrange', tickrange)
                    # Set the formatter
                    # Note some formatters require 'locator' as keyword arg
                    if formatter in ('date', 'concise'):
                        locator = axis.get_major_locator()
                        formatter_kw.setdefault('locator', locator)
                    formatter = axistools.Formatter(
                        formatter, date=date, **formatter_kw)
                    axis.set_major_formatter(formatter)

                # Ensure no out-of-bounds ticks; set_smart_bounds() can fail
                # * Using set_bounds did not work, so instead just turn
                #   locators into fixed version.
                # * Most locators take no arguments in __call__, and some do
                #   not have tick_values method, so we just call them.
                if (bounds is not None
                    or fixticks
                    or isinstance(formatter, mticker.FixedFormatter)
                        or axis.get_scale() == 'cutoff'):
                    if bounds is None:
                        bounds = getattr(self, 'get_' + x + 'lim')()
                    locator = axistools.Locator([
                        x for x in axis.get_major_locator()()
                        if bounds[0] <= x <= bounds[1]
                    ])
                    axis.set_major_locator(locator)
                    locator = axistools.Locator([
                        x for x in axis.get_minor_locator()()
                        if bounds[0] <= x <= bounds[1]
                    ])
                    axis.set_minor_locator(locator)

            # Call parent
            if aspect is not None:
                self.set_aspect(aspect)
            super().format(**kwargs)

    def altx(self, *args, **kwargs):
        # Cannot wrap twiny() because we want to use XYAxes, not
        # matplotlib Axes. Instead use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py  # noqa
        if self._altx_child:
            raise RuntimeError('No more than *two* twin axes!')
        if self._altx_parent:
            raise RuntimeError('This *is* a twin axes!')
        with self.figure._unlock():
            ax = self._make_twin_axes(sharey=self, projection='xy')
        # shared axes must have matching autoscale
        ax.set_autoscaley_on(self.get_autoscaley_on())
        ax.grid(False)
        self._altx_child = ax
        ax._altx_parent = self
        self._altx_overrides()
        ax._altx_overrides()
        self.add_child_axes(ax)
        return ax

    def alty(self):
        if self._alty_child:
            raise RuntimeError('No more than *two* twin axes!')
        if self._alty_parent:
            raise RuntimeError('This *is* a twin axes!')
        with self.figure._unlock():
            ax = self._make_twin_axes(sharex=self, projection='xy')
        # shared axes must have matching autoscale
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.grid(False)
        self._alty_child = ax
        ax._alty_parent = self
        self._alty_overrides()
        ax._alty_overrides()
        self.add_child_axes(ax)
        return ax

    def dualx(self, transform, transform_kw=None, **kwargs):
        # The axis scale is used to transform units on the left axis, linearly
        # spaced, to units on the right axis... so the right scale must scale
        # its data with the *inverse* of this transform. We do this below.
        # NOTE: Matplotlib 3.1 has a 'secondary axis' feature. This one is
        # simpler, because it does not implement the function transform as
        # an axis scale (meaning user just has to supply the forward
        # transformation, not the backwards one), and does not invent a new
        # class with a bunch of complicated setters.
        ax = self.altx()
        funcscale_funcs, funcscale_kw, kwargs = _parse_dualxy_args(
            'x', transform, transform_kw, kwargs
        )
        self._dualx_data = (funcscale_funcs, funcscale_kw)
        self._dualx_overrides()
        ax.format(**kwargs)
        return ax

    def dualy(self, transform, transform_kw=None, **kwargs):
        ax = self.alty()
        funcscale_funcs, funcscale_kw, kwargs = _parse_dualxy_args(
            'y', transform, transform_kw, kwargs
        )
        self._dualy_data = (funcscale_funcs, funcscale_kw)
        self._dualy_overrides()
        ax.format(**kwargs)
        return ax

    def draw(self, renderer=None, *args, **kwargs):
        """Adds post-processing steps before axes is drawn."""
        # NOTE: This mimics matplotlib API, which calls identical
        # post-processing steps in both draw() and get_tightbbox()
        self._hide_labels()
        self._altx_overrides()
        self._alty_overrides()
        self._dualx_overrides()
        self._dualy_overrides()
        self._datex_rotate()
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        super().draw(renderer, *args, **kwargs)

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Adds post-processing steps before tight bounding box is
        calculated."""
        self._hide_labels()
        self._altx_overrides()
        self._alty_overrides()
        self._dualx_overrides()
        self._dualy_overrides()
        self._datex_rotate()
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        return super().get_tightbbox(renderer, *args, **kwargs)

    def twinx(self):
        return self.alty()

    def twiny(self):
        return self.altx()

    altx.__doc__ = altxy_descrip % {
        'x': 'x', 'x1': 'bottom', 'x2': 'top',
        'y': 'y', 'y1': 'left', 'y2': 'right',
    }
    alty.__doc__ = altxy_descrip % {
        'x': 'y', 'x1': 'left', 'x2': 'right',
        'y': 'x', 'y1': 'bottom', 'y2': 'top',
    }
    dualx.__doc__ = dualxy_descrip % {
        'x': 'x', 'args': ', '.join(dualxy_kwargs)
    }
    dualy.__doc__ = dualxy_descrip % {
        'x': 'y', 'args': ', '.join(dualxy_kwargs)
    }
    twinx.__doc__ = twinxy_descrip % {
        'x': 'y', 'x1': 'left', 'x2': 'right',
        'y': 'x', 'y1': 'bottom', 'y2': 'top',
    }
    twiny.__doc__ = twinxy_descrip % {
        'x': 'x', 'x1': 'bottom', 'x2': 'top',
        'y': 'y', 'y1': 'left', 'y2': 'right',
    }


class PolarAxes(Axes, mproj.PolarAxes):
    """Intermediate class, mixes `ProjAxes` with
    `~matplotlib.projections.polar.PolarAxes`."""
    #: The registered projection name.
    name = 'polar'

    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`, `Axes`
        """
        # Set tick length to zero so azimuthal labels are not too offset
        # Change default radial axis formatter but keep default theta one
        super().__init__(*args, **kwargs)
        formatter = axistools.Formatter('auto')
        self.yaxis.set_major_formatter(formatter)
        self.yaxis.isDefault_majfmt = True
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)

    def format(self, *args,
               r0=None, theta0=None, thetadir=None,
               thetamin=None, thetamax=None, thetalim=None,
               rmin=None, rmax=None, rlim=None,
               rlabelpos=None, rscale=None, rborder=None,
               thetalocator=None, rlocator=None, thetalines=None, rlines=None,
               thetaformatter=None, rformatter=None,
               thetalabels=None, rlabels=None,
               thetalocator_kw=None, rlocator_kw=None,
               thetaformatter_kw=None, rformatter_kw=None,
               **kwargs):
        """
        Calls `Axes.format` and `Axes.context`, formats radial gridline
        locations, gridline labels, limits, and more. All ``theta`` arguments
        are specified in *degrees*, not radians. The below parameters are
        specific to `PolarAxes`.

        Parameters
        ----------
        r0 : float, optional
            The radial origin.
        theta0 : {'N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'}
            The zero azimuth location.
        thetadir : {-1, 1, 'clockwise', 'anticlockwise', 'counterclockwise'}, \
optional
            The positive azimuth direction. Clockwise corresponds to ``-1``
            and anticlockwise corresponds to ``-1``. Default is ``-1``.
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

        See also
        --------
        :py:obj:`Axes.format`, :py:obj:`Axes.context`
        """
        context, kwargs = self.context(**kwargs)
        with context:
            # Not mutable default args
            thetalocator_kw = thetalocator_kw or {}
            thetaformatter_kw = thetaformatter_kw or {}
            rlocator_kw = rlocator_kw or {}
            rformatter_kw = rformatter_kw or {}
            # Flexible input
            if rlim is not None:
                if rmin is not None or rmax is not None:
                    warnings.warn(
                        f'Conflicting keyword args rmin={rmin}, rmax={rmax}, '
                        f'and rlim={rlim}. Using "rlim".')
                rmin, rmax = rlim
            if thetalim is not None:
                if thetamin is not None or thetamax is not None:
                    warnings.warn(
                        f'Conflicting keyword args thetamin={thetamin}, '
                        f'thetamax={thetamax}, and thetalim={thetalim}. '
                        f'Using "thetalim".')
                thetamin, thetamax = thetalim
            thetalocator = _notNone(
                thetalines, thetalocator, None,
                names=('thetalines', 'thetalocator'))
            thetaformatter = _notNone(
                thetalabels, thetaformatter, None,
                names=('thetalabels', 'thetaformatter'))
            rlocator = _notNone(rlines, rlocator, None,
                                names=('rlines', 'rlocator'))
            rformatter = _notNone(rlabels, rformatter,
                                  None, names=('rlabels', 'rformatter'))

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
            for (
                x, r, axis,
                min_, max_,
                locator, formatter,
                locator_kw, formatter_kw,
            ) in zip(
                ('x', 'y'), ('theta', 'r'), (self.xaxis, self.yaxis),
                (thetamin, rmin), (thetamax, rmax),
                (thetalocator, rlocator), (thetaformatter, rformatter),
                (thetalocator_kw, rlocator_kw),
                (thetaformatter_kw, rformatter_kw)
            ):
                # Axis limits
                # Try to use public API where possible
                if min_ is not None:
                    getattr(self, 'set_' + r + 'min')(min_)
                else:
                    min_ = getattr(self, 'get_' + r + 'min')()
                if max_ is not None:
                    getattr(self, 'set_' + r + 'max')(max_)
                else:
                    max_ = getattr(self, 'get_' + r + 'max')()

                # Spine settings
                kw = rc.fill({
                    'linewidth': 'axes.linewidth',
                    'color': 'axes.edgecolor',
                })
                sides = ('inner', 'polar') if r == 'r' else ('start', 'end')
                spines = [self.spines[s] for s in sides]
                for spine, side in zip(spines, sides):
                    spine.update(kw)

                # Grid and grid label settings
                # NOTE: Not sure if polar lines inherit tick or grid props
                kw = rc.fill({
                    'color': x + 'tick.color',
                    'labelcolor': 'tick.labelcolor',  # new props
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
                    'weight': 'tick.labelweight'
                })
                for t in axis.get_ticklabels():
                    t.update(kw)

                # Tick locator, which in this case applies to gridlines
                # NOTE: Must convert theta locator input to radians, then back
                # to degrees.
                if locator is not None:
                    if r == 'theta' and (
                            not isinstance(locator, (str, mticker.Locator))):
                        # real axis limts are rad
                        locator = np.deg2rad(locator)
                    locator = axistools.Locator(locator, **locator_kw)
                    locator.set_axis(axis)  # this is what set_locator does
                    grids = np.array(locator())
                    if r == 'r':
                        grids = grids[(grids >= min_) & (grids <= max_)]
                        self.set_rgrids(grids)
                    else:
                        grids = np.rad2deg(grids)
                        grids = grids[(grids >= min_) & (grids <= max_)]
                        if grids[-1] == min_ + 360:  # exclusive if 360 degrees
                            grids = grids[:-1]
                        self.set_thetagrids(grids)
                # Tick formatter and toggling
                if formatter is not None:
                    formatter = axistools.Formatter(formatter, **formatter_kw)
                    axis.set_major_formatter(formatter)

            # Parent method
            super().format(*args, **kwargs)

    # Disabled methods suitable only for cartesian axes
    _disable = _disable_decorator(
        'Invalid plotting method {!r} for polar axes.')
    twinx = _disable(Axes.twinx)
    twiny = _disable(Axes.twiny)
    matshow = _disable(Axes.matshow)
    imshow = _disable(Axes.imshow)
    spy = _disable(Axes.spy)
    hist = _disable(Axes.hist)
    hist2d = _disable(Axes.hist2d)
    boxplot = _disable(Axes.boxplot)
    violinplot = _disable(Axes.violinplot)
    step = _disable(Axes.step)
    stem = _disable(Axes.stem)
    stackplot = _disable(Axes.stackplot)
    table = _disable(Axes.table)
    eventplot = _disable(Axes.eventplot)
    pie = _disable(Axes.pie)
    xcorr = _disable(Axes.xcorr)
    acorr = _disable(Axes.acorr)
    psd = _disable(Axes.psd)
    csd = _disable(Axes.csd)
    cohere = _disable(Axes.cohere)
    specgram = _disable(Axes.specgram)
    angle_spectrum = _disable(Axes.angle_spectrum)
    phase_spectrum = _disable(Axes.phase_spectrum)
    magnitude_spectrum = _disable(Axes.magnitude_spectrum)


def _circle_path(N=100):
    """Return a circle `~matplotlib.path.Path` used as the outline
    for polar stereographic, azimuthal equidistant, and Lambert
    conformal projections. This was developed from `this cartopy example \
<https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`__."""  # noqa
    theta = np.linspace(0, 2 * np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)


class ProjAxes(Axes):
    """Intermediate class, shared by `GeoAxes` and
    `BasemapAxes`. Disables methods that are inappropriate for map
    projections and adds `ProjAxes.format`, so that arguments
    passed to `Axes.format` are identical for `GeoAxes`
    and `BasemapAxes`."""

    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`, `Axes`, `GeoAxes`, `BasemapAxes`
        """
        # Store props that let us dynamically and incrementally modify
        # line locations and settings like with Cartesian axes
        self._boundinglat = None
        self._latmax = None
        self._latlines = None
        self._lonlines = None
        self._lonlines_values = None
        self._latlines_values = None
        self._lonlines_labels = None
        self._latlines_labels = None
        super().__init__(*args, **kwargs)

    def format(self, *,
               lonlim=None, latlim=None, boundinglat=None, grid=None,
               lonlines=None, lonlocator=None,
               latlines=None, latlocator=None, latmax=None,
               labels=None, latlabels=None, lonlabels=None,
               patch_kw=None, **kwargs,
               ):
        """
        Calls `Axes.format` and `Axes.context`, formats the meridian
        and parallel labels, longitude and latitude map limits, geographic
        features, and more.

        Parameters
        ----------
        lonlim, latlim : (float, float), optional
            Longitude and latitude limits of projection, applied
            with `~cartopy.mpl.geoaxes.GeoAxes.set_extent`.
            For cartopy axes only.
        boundinglat : float, optional
            The edge latitude for the circle bounding North Pole and
            South Pole-centered projections. For cartopy axes only.
        grid : bool, optional
            Toggles meridian and parallel gridlines on and off. Default is
            :rc:`geogrid`.
        lonlines, latlines : float or list of float, optional
            If float, indicates the *spacing* of meridian and parallel
            gridlines. Otherwise, must be a list of floats indicating specific
            meridian and parallel gridlines to draw.
        lonlocator, latlocator : optional
            Aliases for `lonlines`, `latlines`.
        latmax : float, optional
            The maximum absolute latitude for meridian gridlines. Default is
            :rc:`geogrid.latmax`.
        labels : bool, optional
            Toggles meridian and parallel gridline labels on and off. Default
            is :rc:`geogrid.labels`.
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

        land, ocean, coast, rivers, lakes, borders, innerborders : bool, \
optional
            Toggles various geographic features. These are actually the
            :rcraw:`land`, :rcraw:`ocean`, :rcraw:`coast`, :rcraw:`rivers`,
            :rcraw:`lakes`, :rcraw:`borders`, and :rcraw:`innerborders`
            settings passed to `~proplot.axes.Axes.context`. The style can
            be modified by passing additional settings, e.g.
            :rcraw:`landcolor`.
        patch_kw : dict-like, optional
            Keyword arguments used to update the background patch object. You
            can use this, for example, to set background hatching with
            ``patch_kw={'hatch':'xxx'}``.
        **kwargs
            Passed to `Axes.format` and `Axes.context`.

        See also
        --------
        :py:obj:`Axes.format`, :py:obj:`Axes.context`
        """
        # Parse alternative keyword args
        # TODO: Why isn't default latmax 80 respected sometimes?
        context, kwargs = self.context(**kwargs)
        with context:
            lonlines = _notNone(
                lonlines, lonlocator, rc['geogrid.lonstep'],
                names=('lonlines', 'lonlocator'))
            latlines = _notNone(
                latlines, latlocator, rc['geogrid.latstep'],
                names=('latlines', 'latlocator'))
            latmax = _notNone(latmax, rc['geogrid.latmax'])
            labels = _notNone(labels, rc['geogrid.labels'])
            grid = _notNone(grid, rc['geogrid'])
            if labels:
                lonlabels = _notNone(lonlabels, 1)
                latlabels = _notNone(latlabels, 1)

            # Longitude gridlines, draw relative to projection prime meridian
            # NOTE: Always generate gridlines array on first format call
            # because rc setting will be not None
            if isinstance(self, GeoAxes):
                lon_0 = self.projection.proj4_params.get('lon_0', 0)
            else:
                base = 5
                lon_0 = base * round(
                    self.projection.lonmin / base) + 180  # central longitude
            if lonlines is not None:
                if not np.iterable(lonlines):
                    lonlines = utils.arange(lon_0 - 180, lon_0 + 180, lonlines)
                    lonlines = lonlines.astype(np.float64)
                    lonlines[-1] -= 1e-10  # make sure appears on *right*
                lonlines = [*lonlines]

            # Latitudes gridlines, draw from -latmax to latmax unless result
            # would be asymmetrical across equator
            # NOTE: Basemap axes redraw *meridians* if they detect latmax was
            # explicitly changed, so important not to overwrite 'latmax'
            # with default value! Just need it for this calculation, then when
            # drawparallels is called will use self._latmax
            if latlines is not None or latmax is not None:
                # Fill defaults
                if latlines is None:
                    latlines = _notNone(
                        self._latlines_values, rc.get('geogrid.latstep'))
                ilatmax = _notNone(latmax, self._latmax,
                                   rc.get('geogrid.latmax'))
                # Get tick locations
                if not np.iterable(latlines):
                    if (ilatmax % latlines) == (-ilatmax % latlines):
                        latlines = utils.arange(-ilatmax, ilatmax, latlines)
                    else:
                        latlines = utils.arange(0, ilatmax, latlines)
                        if latlines[-1] != ilatmax:
                            latlines = np.concatenate((latlines, [ilatmax]))
                        latlines = np.concatenate(
                            (-latlines[::-1], latlines[1:]))
                latlines = [*latlines]

            # Length-4 boolean arrays of whether and where to toggle labels
            # Format is [left, right, bottom, top]
            lonarray, latarray = [], []
            for labs, array in zip(
                (lonlabels, latlabels), (lonarray, latarray)
            ):
                if labs is None:
                    continue  # leave empty
                if isinstance(labs, str):
                    string = labs
                    labs = [0] * 4
                    for idx, char in zip([0, 1, 2, 3], 'lrbt'):
                        if char in string:
                            labs[idx] = 1
                elif not np.iterable(labs):
                    labs = np.atleast_1d(labs)
                if len(labs) == 1:
                    labs = [*labs, 0]  # default is to label bottom/left
                if len(labs) == 2:
                    if array is lonarray:
                        labs = [0, 0, *labs]
                    else:
                        labs = [*labs, 0, 0]
                elif len(labs) != 4:
                    raise ValueError(f'Invalid lon/lat label spec: {labs}.')
                array[:] = labs
            lonarray = lonarray or None  # None so use default locations
            latarray = latarray or None

            # Add attributes for redrawing lines
            if latmax is not None:
                self._latmax = latmax
            if latlines is not None:
                self._latlines_values = latlines
            if lonlines is not None:
                self._lonlines_values = lonlines
            if latarray is not None:
                self._latlines_labels = latarray
            if lonarray is not None:
                self._lonlines_labels = lonarray

            # Grid toggling, must come after everything else in case e.g.
            # rc.geogrid is False but user passed grid=True so we need to
            # recover the *default* lonlines and latlines values
            if grid is not None:
                if not grid:
                    lonlines = latlines = []
                else:
                    lonlines = self._lonlines_values
                    latlines = self._latlines_values

            # Apply formatting to basemap or cartpoy axes
            patch_kw = patch_kw or {}
            self._format_apply(patch_kw, lonlim, latlim, boundinglat,
                               lonlines, latlines, latmax, lonarray, latarray)
            super().format(**kwargs)

    # Disabled methods suitable only for cartesian axes
    _disable = _disable_decorator(
        'Invalid plotting method {!r} for map projection axes.')
    bar = _disable(Axes.bar)
    barh = _disable(Axes.barh)
    twinx = _disable(Axes.twinx)
    twiny = _disable(Axes.twiny)
    matshow = _disable(Axes.matshow)
    imshow = _disable(Axes.imshow)
    spy = _disable(Axes.spy)
    hist = _disable(Axes.hist)
    hist2d = _disable(Axes.hist2d)
    boxplot = _disable(Axes.boxplot)
    violinplot = _disable(Axes.violinplot)
    step = _disable(Axes.step)
    stem = _disable(Axes.stem)
    stackplot = _disable(Axes.stackplot)
    table = _disable(Axes.table)
    eventplot = _disable(Axes.eventplot)
    pie = _disable(Axes.pie)
    xcorr = _disable(Axes.xcorr)
    acorr = _disable(Axes.acorr)
    psd = _disable(Axes.psd)
    csd = _disable(Axes.csd)
    cohere = _disable(Axes.cohere)
    specgram = _disable(Axes.specgram)
    angle_spectrum = _disable(Axes.angle_spectrum)
    phase_spectrum = _disable(Axes.phase_spectrum)
    magnitude_spectrum = _disable(Axes.magnitude_spectrum)


def _add_gridline_label(self, value, axis, upper_end):
    """Gridliner method monkey patch. Always print number in range
    (180W, 180E)."""
    # Have 3 choices (see Issue #78):
    # 1. lonlines go from -180 to 180, but get double 180 labels at dateline
    # 2. lonlines go from -180 to e.g. 150, but no lines from 150 to dateline
    # 3. lonlines go from lon_0 - 180 to lon_0 + 180 mod 360, but results
    #    in non-monotonic array causing double gridlines east of dateline
    # 4. lonlines go from lon_0 - 180 to lon_0 + 180 monotonic, but prevents
    #    labels from being drawn outside of range (-180, 180)
    # These monkey patches choose #4 and permit labels being drawn
    # outside of (-180 180)
    if axis == 'x':
        value = (value + 180) % 360 - 180
    return type(self)._add_gridline_label(self, value, axis, upper_end)


def _axes_domain(self, *args, **kwargs):
    """Gridliner method monkey patch. Filter valid label coordinates to values
    between lon_0 - 180 and lon_0 + 180."""
    # See _add_gridline_label for detials
    lon_0 = self.axes.projection.proj4_params.get('lon_0', 0)
    x_range, y_range = type(self)._axes_domain(self, *args, **kwargs)
    x_range = np.asarray(x_range) + lon_0
    return x_range, y_range


class GeoAxes(ProjAxes, GeoAxes):
    """Axes subclass for plotting `cartopy \
<https://scitools.org.uk/cartopy/docs/latest/>`__ projections. Initializes
    the `cartopy.crs.Projection` instance, enforces `global extent \
<https://stackoverflow.com/a/48956844/4970632>`__
    for most projections by default, and draws `circular boundaries \
<https://scitools.org.uk/cartopy/docs/latest/gallery/always_circular_stereo.html>`__
    around polar azimuthal, stereographic, and Gnomonic projections bounded at
    the equator by default."""  # noqa
    #: The registered projection name.
    name = 'geo'
    _n_points = 100  # number of points for drawing circle map boundary

    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~cartopy.crs.Projection`
            The `~cartopy.crs.Projection` instance.
        *args, **kwargs
            Passed to `~cartopy.mpl.geoaxes.GeoAxes`.

        See also
        --------
        `~proplot.subplots.subplots`, `Axes`, `~proplot.projs.Proj`
        """
        # GeoAxes initialization. Note that critical attributes like
        # outline_patch needed by _format_apply are added before it is called.
        import cartopy.crs as ccrs
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError(
                'GeoAxes requires map_projection=cartopy.crs.Projection.')
        super().__init__(*args, map_projection=map_projection, **kwargs)

        # Zero out ticks so gridlines are not offset
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)

        # Set extent and boundary extent for projections
        # The default bounding latitude is set in _format_apply
        # NOTE: set_global does not mess up non-global projections like OSNI
        if isinstance(self.projection, (
                ccrs.NorthPolarStereo, ccrs.SouthPolarStereo,
                projs.NorthPolarGnomonic, projs.SouthPolarGnomonic,
                projs.NorthPolarAzimuthalEquidistant,
                projs.NorthPolarLambertAzimuthalEqualArea,
                projs.SouthPolarAzimuthalEquidistant,
                projs.SouthPolarLambertAzimuthalEqualArea)):
            self.set_boundary(_circle_path(100), transform=self.transAxes)
        else:
            self.set_global()

    def _format_apply(self, patch_kw, lonlim, latlim, boundinglat,
                      lonlines, latlines, latmax, lonarray, latarray):
        """Apply formatting to cartopy axes."""
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs
        from cartopy.mpl import gridliner

        # Initial gridliner object, which ProPlot passively modifies
        # TODO: Flexible formatter?
        if not self._gridliners:
            gl = self.gridlines(zorder=2.5)  # below text only
            gl._axes_domain = _axes_domain.__get__(gl)  # apply monkey patches
            gl._add_gridline_label = _add_gridline_label.__get__(gl)
            gl.xlines = False
            gl.ylines = False
            try:
                lonformat = gridliner.LongitudeFormatter  # newer
                latformat = gridliner.LatitudeFormatter
            except AttributeError:
                lonformat = gridliner.LONGITUDE_FORMATTER  # older
                latformat = gridliner.LATITUDE_FORMATTER
            gl.xformatter = lonformat
            gl.yformatter = latformat
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            gl.ylabels_left = False
            gl.ylabels_right = False

        # Projection extent
        # NOTE: They may add this as part of set_xlim and set_ylim in future
        # See: https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638  # noqa
        # WARNING: The set_extent method tries to set a *rectangle* between
        # the *4* (x,y) coordinate pairs (each corner), so something like
        # (-180,180,-90,90) will result in *line*, causing error!
        proj = self.projection.proj4_params['proj']
        north = isinstance(self.projection, (
            ccrs.NorthPolarStereo, projs.NorthPolarGnomonic,
            projs.NorthPolarAzimuthalEquidistant,
            projs.NorthPolarLambertAzimuthalEqualArea))
        south = isinstance(self.projection, (
            ccrs.SouthPolarStereo, projs.SouthPolarGnomonic,
            projs.SouthPolarAzimuthalEquidistant,
            projs.SouthPolarLambertAzimuthalEqualArea))
        if north or south:
            if (lonlim is not None or latlim is not None):
                warnings.warn(
                    f'{proj!r} extent is controlled by "boundinglat", '
                    f'ignoring lonlim={lonlim!r} and latlim={latlim!r}.')
            if self._boundinglat is None:
                if isinstance(self.projection, projs.NorthPolarGnomonic):
                    boundinglat = 30
                elif isinstance(self.projection, projs.SouthPolarGnomonic):
                    boundinglat = -30
                else:
                    boundinglat = 0
            if boundinglat is not None and boundinglat != self._boundinglat:
                eps = 1e-10  # bug with full -180, 180 range when lon_0 != 0
                lat0 = (90 if north else -90)
                lon0 = self.projection.proj4_params.get('lon_0', 0)
                extent = [lon0 - 180 + eps,
                          lon0 + 180 - eps, boundinglat, lat0]
                self.set_extent(extent, crs=ccrs.PlateCarree())
                self._boundinglat = boundinglat
        else:
            if boundinglat is not None:
                warnings.warn(
                    f'{proj!r} extent is controlled by "lonlim" and "latlim", '
                    f'ignoring boundinglat={boundinglat!r}.')
            if lonlim is not None or latlim is not None:
                lonlim = lonlim or [None, None]
                latlim = latlim or [None, None]
                lonlim, latlim = [*lonlim], [*latlim]
                lon_0 = self.projection.proj4_params.get('lon_0', 0)
                if lonlim[0] is None:
                    lonlim[0] = lon_0 - 180
                if lonlim[1] is None:
                    lonlim[1] = lon_0 + 180
                eps = 1e-10  # bug with full -180, 180 range when lon_0 != 0
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
        # Collection props, see GoeAxes.gridlines() source code
        kw = rc.fill({
            'alpha': 'geogrid.alpha',
            'color': 'geogrid.color',
            'linewidth': 'geogrid.linewidth',
            'linestyle': 'geogrid.linestyle',
        })  # cached changes
        gl.collection_kwargs.update(kw)
        # Grid locations
        eps = 1e-10
        if lonlines is not None:
            if len(lonlines) == 0:
                gl.xlines = False
            else:
                gl.xlines = True
                gl.xlocator = mticker.FixedLocator(lonlines)
        if latlines is not None:
            if len(latlines) == 0:
                gl.ylines = False
            else:
                gl.ylines = True
                if latlines[0] == -90:
                    latlines[0] += eps
                if latlines[-1] == 90:
                    latlines[-1] -= eps
                gl.ylocator = mticker.FixedLocator(latlines)
        # Grid label toggling
        # Issue warning instead of error!
        if not isinstance(self.projection, (ccrs.Mercator, ccrs.PlateCarree)):
            if latarray is not None and any(latarray):
                warnings.warn(
                    'Cannot add gridline labels to cartopy '
                    f'{type(self.projection).__name__} projection.')
                latarray = [0] * 4
            if lonarray is not None and any(lonarray):
                warnings.warn(
                    'Cannot add gridline labels to cartopy '
                    f'{type(self.projection).__name__} projection.')
                lonarray = [0] * 4
        if latarray is not None:
            gl.ylabels_left = latarray[0]
            gl.ylabels_right = latarray[1]
        if lonarray is not None:
            gl.xlabels_bottom = lonarray[2]
            gl.xlabels_top = lonarray[3]

        # Geographic features
        # WARNING: Seems cartopy features can't be updated!
        # See: https://scitools.org.uk/cartopy/docs/v0.14/_modules/cartopy/feature.html#Feature  # noqa
        # Change the _kwargs property also does *nothing*
        # WARNING: Changing linewidth is impossible with cfeature. Bug?
        # See: https://stackoverflow.com/questions/43671240/changing-line-width-of-cartopy-borders  # noqa
        # TODO: Editing existing natural features? Creating natural features
        # at __init__ time and hiding them?
        # NOTE: The natural_earth_shp method is deprecated, use add_feature.
        # See: https://cartopy-pelson.readthedocs.io/en/readthedocs/whats_new.html  # noqa
        # NOTE: The e.g. cfeature.COASTLINE features are just for convenience,
        # hi res versions. Use cfeature.COASTLINE.name to see how it can be
        # looked up with NaturalEarthFeature.
        reso = rc.get('reso')
        if reso not in ('lo', 'med', 'hi'):
            raise ValueError(f'Invalid resolution {reso}.')
        reso = {
            'lo': '110m',
            'med': '50m',
            'hi': '10m',
        }.get(reso)
        features = {
            'land': ('physical', 'land'),
            'ocean': ('physical', 'ocean'),
            'lakes': ('physical', 'lakes'),
            'coast': ('physical', 'coastline'),
            'rivers': ('physical', 'rivers_lake_centerlines'),
            'borders': ('cultural', 'admin_0_boundary_lines_land'),
            'innerborders': ('cultural', 'admin_1_states_provinces_lakes'),
        }
        for name, args in features.items():
            # Get feature
            if not rc.get(name):  # toggled
                continue
            if getattr(self, '_' + name, None):  # already drawn
                continue
            feat = cfeature.NaturalEarthFeature(*args, reso)
            # For 'lines', need to specify edgecolor and facecolor
            # See: https://github.com/SciTools/cartopy/issues/803
            kw = rc.category(name, cache=False)
            if name in ('coast', 'rivers', 'borders', 'innerborders'):
                kw['edgecolor'] = kw.pop('color')
                kw['facecolor'] = 'none'
            else:
                kw['linewidth'] = 0
            if name in ('ocean',):
                kw['zorder'] = 0.5  # below everything!
            self.add_feature(feat, **kw)
            setattr(self, '_' + name, feat)

        # Update patch
        kw_face = rc.fill({
            'facecolor': 'geoaxes.facecolor'
        })
        kw_face.update(patch_kw)
        self.background_patch.update(kw_face)
        kw_edge = rc.fill({
            'edgecolor': 'geoaxes.edgecolor',
            'linewidth': 'geoaxes.linewidth'
        })
        self.outline_patch.update(kw_edge)

    def _hide_labels(self):
        """No-op for now. In future will hide meridian and parallel labels
        for rectangular projections."""
        pass

    def get_tightbbox(self, renderer, *args, **kwargs):
        """Draw gridliner objects so tight bounding box algorithm will
        incorporate gridliner labels."""
        self._hide_labels()
        if self.get_autoscale_on() and self.ignore_existing_data_limits:
            self.autoscale_view()
        if self.background_patch.reclip:
            clipped_path = self.background_patch.orig_path.clip_to_bbox(
                self.viewLim)
            self.background_patch._path = clipped_path
        self.apply_aspect()
        for gl in self._gridliners:
            try:  # new versions only
                gl._draw_gridliner(background_patch=self.background_patch,
                                   renderer=renderer)
            except TypeError:
                gl._draw_gridliner(background_patch=self.background_patch)
        self._gridliners = []
        return super().get_tightbbox(renderer, *args, **kwargs)

    # Document projection property
    @property
    def projection(self):
        """The `~cartopy.crs.Projection` instance associated with this axes."""
        return self._map_projection

    @projection.setter
    def projection(self, map_projection):
        self._map_projection = map_projection

    # Wrapped methods
    # TODO: Remove this duplication of Axes! Can do this when we implement
    # all wrappers as decorators.
    if GeoAxes is not object:
        text = _text_wrapper(
            GeoAxes.text
        )
        # Wrapped by standardize method
        plot = _default_transform(_plot_wrapper(_standardize_1d(
            _add_errorbars(_cycle_changer(GeoAxes.plot))
        )))
        scatter = _default_transform(_scatter_wrapper(_standardize_1d(
            _add_errorbars(_cycle_changer(GeoAxes.scatter))
        )))
        fill_between = _fill_between_wrapper(_standardize_1d(_cycle_changer(
            GeoAxes.fill_between
        )))
        fill_betweenx = _fill_betweenx_wrapper(_standardize_1d(_cycle_changer(
            GeoAxes.fill_betweenx
        )))
        contour = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.contour
        )))
        contourf = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.contourf
        )))
        pcolor = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.pcolor
        )))
        pcolormesh = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.pcolormesh
        )))
        quiver = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.quiver
        )))
        streamplot = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.streamplot
        )))
        barbs = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxes.barbs
        )))

        # Wrapped only by cmap wrapper
        tripcolor = _default_transform(_cmap_changer(
            GeoAxes.tripcolor
        ))
        tricontour = _default_transform(_cmap_changer(
            GeoAxes.tricontour
        ))
        tricontourf = _default_transform(_cmap_changer(
            GeoAxes.tricontourf
        ))

        # Special GeoAxes commands
        get_extent = _default_crs(GeoAxes.get_extent)
        set_extent = _default_crs(GeoAxes.set_extent)
        set_xticks = _default_crs(GeoAxes.set_xticks)
        set_yticks = _default_crs(GeoAxes.set_yticks)


class BasemapAxes(ProjAxes):
    """Axes subclass for plotting `~mpl_toolkits.basemap` projections. The
    `~mpl_toolkits.basemap.Basemap` projection instance is added as
    the `map_projection` attribute, but this is all abstracted away. You can
    use `~matplotlib.axes.Axes` methods like `~matplotlib.axes.Axes.plot` and
    `~matplotlib.axes.Axes.contour` with your raw longitude-latitude data."""
    #: The registered projection name.
    name = 'basemap'
    _proj_non_rectangular = (
        'ortho', 'geos', 'nsper',
        'moll', 'hammer', 'robin',
        'eck4', 'kav7', 'mbtfpq',
        'sinu', 'vandg',
        'npstere', 'spstere', 'nplaea',
        'splaea', 'npaeqd', 'spaeqd',
    )  # do not use axes spines as boundaries

    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        **kwargs
            Passed to `Axes`.

        See also
        --------
        `~proplot.subplots.subplots`, `Axes`, `~proplot.projs.Proj`
        """
        # Map boundary notes
        # * Must set boundary before-hand, otherwise the set_axes_limits method
        #   called by mcontourf/mpcolormesh/etc draws two mapboundary Patch
        #   objects called "limb1" and "limb2" automatically: one for fill and
        #   the other for the edges
        # * Then, since the patch object in _mapboundarydrawn is only the
        #   fill-version, calling drawmapboundary again will replace only *that
        #   one*, but the original visible edges are still drawn -- so e.g. you
        #   can't change the color
        # * If you instead call drawmapboundary right away, _mapboundarydrawn
        #   will contain both the edges and the fill; so calling it again will
        #   replace *both*
        import mpl_toolkits.basemap as mbasemap  # verify package is available
        if not isinstance(map_projection, mbasemap.Basemap):
            raise ValueError(
                'BasemapAxes requires map_projection=basemap.Basemap')
        self._map_projection = map_projection
        self._map_boundary = None
        self._has_recurred = False  # use this to override plotting methods
        super().__init__(*args, **kwargs)

    def _format_apply(self, patch_kw, lonlim, latlim, boundinglat,
                      lonlines, latlines, latmax, lonarray, latarray):
        """Applies formatting to basemap axes."""
        # Checks
        if (lonlim is not None or latlim is not None
                or boundinglat is not None):
            warnings.warn(f'Got lonlim={lonlim!r}, latlim={latlim!r}, '
                          f'boundinglat={boundinglat!r}, but you cannot "zoom '
                          'into" a basemap projection after creating it. '
                          'Pass proj_kw in your call to subplots '
                          'with any of the following basemap keywords: '
                          "'boundinglat', 'llcrnrlon', 'llcrnrlat', "
                          "'urcrnrlon', 'urcrnrlat', 'llcrnrx', 'llcrnry', "
                          "'urcrnrx', 'urcrnry', 'width', or 'height'.")

        # Map boundary
        # * First have to *manually replace* the old boundary by just
        #   deleting the original one
        # * If boundary is drawn successfully should be able to call
        #   self.projection._mapboundarydrawn.set_visible(False) and
        #   edges/fill color disappear
        # * For now will enforce that map plots *always* have background
        #   whereas axes plots can have transparent background
        kw_edge = rc.fill({
            'linewidth': 'geoaxes.linewidth',
            'edgecolor': 'geoaxes.edgecolor'
        })
        kw_face = rc.fill({
            'facecolor': 'geoaxes.facecolor'
        })
        patch_kw = patch_kw or {}
        kw_face.update(patch_kw)
        self.axesPatch = self.patch  # bugfix or something
        if self.projection.projection in self._proj_non_rectangular:
            self.patch.set_alpha(0)  # make patch invisible
            if not self.projection._mapboundarydrawn:
                # set fill_color to 'none' to make transparent
                p = self.projection.drawmapboundary(ax=self)
            else:
                p = self.projection._mapboundarydrawn
            p.update({**kw_face, **kw_edge})
            p.set_rasterized(False)
            p.set_clip_on(False)  # so edges denoting boundary aren't cut off
            self._map_boundary = p
        else:
            self.patch.update({**kw_face, 'edgecolor': 'none'})
            for spine in self.spines.values():
                spine.update(kw_edge)

        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian
        # markers come up to the axes bbox
        lkw = rc.fill({
            'alpha': 'geogrid.alpha',
            'color': 'geogrid.color',
            'linewidth': 'geogrid.linewidth',
            'linestyle': 'geogrid.linestyle',
        }, cache=False)
        tkw = rc.fill({
            'color': 'geogrid.color',
            'fontsize': 'geogrid.labelsize',
        }, cache=False)
        # Change from left/right/bottom/top to left/right/top/bottom
        if lonarray is not None:
            lonarray[2:] = lonarray[2:][::-1]
        if latarray is not None:
            latarray[2:] = latarray[2:][::-1]

        # Parallel lines
        if latlines is not None or latmax is not None or latarray is not None:
            if self._latlines:
                for pi in self._latlines.values():
                    for obj in [i for j in pi for i in j]:  # magic
                        obj.set_visible(False)
            ilatmax = _notNone(latmax, self._latmax)
            latlines = _notNone(latlines, self._latlines_values)
            latarray = _notNone(latarray, self._latlines_labels, [0] * 4)
            p = self.projection.drawparallels(
                latlines, latmax=ilatmax, labels=latarray, ax=self)
            for pi in p.values():  # returns dict, where each one is tuple
                # Tried passing clip_on to the below, but it does nothing
                # Must set for lines created after the fact
                for obj in [i for j in pi for i in j]:
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
            self._latlines = p

        # Meridian lines
        if lonlines is not None or latmax is not None or lonarray is not None:
            if self._lonlines:
                for pi in self._lonlines.values():
                    for obj in [i for j in pi for i in j]:  # magic
                        obj.set_visible(False)
            ilatmax = _notNone(latmax, self._latmax)
            lonlines = _notNone(lonlines, self._lonlines_values)
            lonarray = _notNone(lonarray, self._lonlines_labels, [0] * 4)
            p = self.projection.drawmeridians(
                lonlines, latmax=ilatmax, labels=lonarray, ax=self)
            for pi in p.values():
                for obj in [i for j in pi for i in j]:
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
            self._lonlines = p

        # Geography
        # TODO: Allow setting the zorder.
        # NOTE: Also notable are drawcounties, blumarble, drawlsmask,
        # shadedrelief, and etopo methods.
        features = {
            'land': 'fillcontinents',
            'coast': 'drawcoastlines',
            'rivers': 'drawrivers',
            'borders': 'drawcountries',
            'innerborders': 'drawstates',
        }
        for name, method in features.items():
            if not rc.get(name):  # toggled
                continue
            if getattr(self, f'_{name}', None):  # already drawn
                continue
            kw = rc.category(name, cache=False)
            feat = getattr(self.projection, method)(ax=self)
            if isinstance(feat, (list, tuple)):  # list of artists?
                for obj in feat:
                    obj.update(kw)
            else:
                feat.update(kw)
            setattr(self, '_' + name, feat)

    # Document projection property
    @property
    def projection(self):
        """The `~mpl_toolkits.basemap.Basemap` instance associated with
        this axes."""
        return self._map_projection

    @projection.setter
    def projection(self, map_projection):
        self._map_projection = map_projection

    # Wrapped methods
    plot = _norecurse(_default_latlon(_plot_wrapper(_standardize_1d(
        _add_errorbars(_cycle_changer(_redirect(maxes.Axes.plot)))
    ))))
    scatter = _norecurse(_default_latlon(_scatter_wrapper(_standardize_1d(
        _add_errorbars(_cycle_changer(_redirect(maxes.Axes.scatter)))
    ))))
    contour = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.contour)
    ))))
    contourf = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.contourf)
    ))))
    pcolor = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.pcolor)
    ))))
    pcolormesh = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.pcolormesh)
    ))))
    quiver = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.quiver)
    ))))
    streamplot = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.streamplot)
    ))))
    barbs = _norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _redirect(maxes.Axes.barbs)
    ))))
    hexbin = _norecurse(_standardize_1d(_cmap_changer(
        _redirect(maxes.Axes.hexbin)
    )))
    imshow = _norecurse(_cmap_changer(
        _redirect(maxes.Axes.imshow)
    ))


# Register the projections
mproj.register_projection(PolarAxes)
mproj.register_projection(XYAxes)
mproj.register_projection(GeoAxes)
mproj.register_projection(BasemapAxes)
