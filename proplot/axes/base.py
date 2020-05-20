#!/usr/bin/env python3
"""
The base axes class used for all ProPlot figures.
"""
import numpy as np
import copy
import matplotlib.axes as maxes
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections
import matplotlib.projections as mprojections
from numbers import Integral, Number
from .plot import (
    _get_transform,
    _bar_wrapper, _barh_wrapper, _boxplot_wrapper,
    _cmap_changer, _cycle_changer,
    _fill_between_wrapper, _fill_betweenx_wrapper, _hlines_wrapper,
    _hist_wrapper, _indicate_error,
    _parametric_wrapper, _plot_wrapper, _scatter_wrapper, _stem_wrapper,
    _standardize_1d, _standardize_2d,
    _text_wrapper, _violinplot_wrapper, _vlines_wrapper,
    colorbar_wrapper, legend_wrapper,
)
from .. import gridspec as pgridspec
from .. import constructor
from ..config import rc
from ..utils import units, edges
from ..internals import ic  # noqa: F401
from ..internals import docstring, rcsetup, warnings, _not_none

__all__ = ['Axes']

ABC_STRING = 'abcdefghijklmnopqrstuvwxyz'
LOC_TRANSLATE = {  # for inset colorbars and legends TODO: also as text locations
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


docstring.snippets['axes.other'] = """
rc_kw : dict, optional
    Dictionary containing `~proplot.config.rc` settings applied to
    this axes using `~proplot.config.RcConfigurator.context`.
**kwargs
    Passed to `Axes.format` or passed to `~proplot.config.RcConfigurator.context`
    and used to update axes `~proplot.config.rc` settings. For example,
    ``abcstyle='A.'`` modifies the :rcraw:`abc.style` setting.
"""

docstring.snippets['axes.patch_kw'] = """
patch_kw : dict-like, optional
    Keyword arguments used to update the background patch object. You
    can use this, for example, to set background hatching with
    ``patch_kw={'hatch': 'xxx'}``.
"""

docstring.snippets['axes.proj'] = """
The map projection specification(s). If ``'cartesian'`` (the default), a
`~proplot.axes.CartesianAxes` is created. If ``'polar'``, a
`~proplot.axes.PolarAxes` is created. Otherwise, the argument is
interpreted by `~proplot.constructor.Proj`, and the result is used
to make a `~proplot.axes.GeoAxes` (in this case the argument can be
a `cartopy.crs.Projection` instance, a `~mpl_toolkits.basemap.Basemap`
instance, or a projection name listed in :ref:`this table <proj_table>`).
"""

docstring.snippets['axes.inset'] = """
Return an inset `CartesianAxes`. This is similar to the builtin
`~matplotlib.axes.Axes.inset_axes` but includes some extra options.

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
proj, projection : str, `cartopy.crs.Projection`, or `~mpl_toolkits.basemap.Basemap`
    The map projection specification(s). If not provided, the inset axes
    projection is identical to the current axes projection. If ``'cartesian'``,
    a `~proplot.axes.CartesianAxes` inset is created. If ``'polar'``, a
    `~proplot.axes.PolarAxes` inset is created. Otherwise, the argument is
    interpreted by `~proplot.constructor.Proj`, and the result is used
    to make a `~proplot.axes.GeoAxes` (in this case the argument can be
    a `cartopy.crs.Projection` instance, a `~mpl_toolkits.basemap.Basemap`
    instance, or a projection name listed in :ref:`this table <proj_table>`).
proj_kw, projection_kw : dict-like, optional
    Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
    cartopy `~cartopy.crs.Projection` classes on instantiation.
basemap : bool or dict-like, optional
    Whether to use `~mpl_toolkits.basemap.Basemap` or
    `~cartopy.crs.Projection` for map projections. Default is ``False``.
zorder : float, optional
    The `zorder <https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
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
    Passed to `CartesianAxes`.
""" % docstring.snippets

docstring.snippets['axes.panel'] = """
Return a panel drawn along the edge of this axes.

Parameters
----------
side : str, optional
    The panel location. The following location keys are valid:

    ==========  =====================
    Location    Valid keys
    ==========  =====================
    left        ``'left'``, ``'l'``
    right       ``'right'``, ``'r'``
    bottom      ``'bottom'``, ``'b'``
    top         ``'top'``, ``'t'``
    ==========  =====================

width : float or str or list thereof, optional
    The panel width. Units are interpreted by `~proplot.utils.units`.
    Default is :rc:`subplots.panelwidth`.
space : float or str or list thereof, optional
    Empty space between the main subplot and the panel.
    When :rcraw:`tight` is ``True``, this is adjusted automatically.
    Otherwise, the default is :rc:`subplots.panelpad`.
share : bool, optional
    Whether to enable axis sharing between the *x* and *y* axes of the
    main subplot and the panel long axes for each panel in the stack.
    Sharing between the panel short axis and other panel short axes
    is determined by figure-wide `sharex` and `sharey` settings.

Returns
-------
`~proplot.axes.CartesianAxes`
    The panel axes.
"""


class Axes(maxes.Axes):
    """
    Lowest-level axes subclass. Handles titles and axis
    sharing. Adds several new methods and overrides existing ones.
    """
    def __init__(self, *args, number=None, main=False, **kwargs):
        """
        Parameters
        ----------
        number : int
            The subplot number, used for a-b-c labeling. See `~Axes.format`
            for details. Note the first axes is ``1``, not ``0``.
        main : bool, optional
            Used internally, indicates whether this is a "main axes" rather
            than a twin, panel, or inset axes.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~matplotlib.axes.Axes`.

        See also
        --------
        matplotlib.axes.Axes
        proplot.axes.CartesianAxes
        proplot.axes.PolarAxes
        proplot.axes.GeoAxes
        """
        super().__init__(*args, **kwargs)

        # Ensure isDefault_minloc enabled at start, needed for dual axes
        self.xaxis.isDefault_minloc = self.yaxis.isDefault_minloc = True

        # Properties
        self._abc_loc = None
        self._abc_text = None
        self._abc_border_kwargs = {}  # abs border properties
        self._title_loc = None  # location of main title
        self._title_pad = rc['axes.titlepad']  # format() can overwrite
        self._title_pad_active = None
        self._title_border_kwargs = {}  # title border properties
        self._above_top_panels = True  # TODO: add rc prop?
        self._bottom_panels = []
        self._top_panels = []
        self._left_panels = []
        self._right_panels = []
        self._tightbbox = None  # bounding boxes are saved
        self._panel_hidden = False  # True when "filled" with cbar/legend
        self._panel_parent = None
        self._panel_share = False
        self._panel_sharex_group = False
        self._panel_sharey_group = False
        self._panel_side = None
        self._inset_parent = None
        self._inset_zoom = False
        self._inset_zoom_data = None
        self._alty_child = None
        self._altx_child = None
        self._alty_parent = None
        self._altx_parent = None
        self.number = number  # for abc numbering
        if main:
            self.figure._axes_main.append(self)

        # On-the-fly legends and colorbars
        self._auto_colorbar = {}
        self._auto_legend = {}

        # Figure row and column labels
        # NOTE: Most of these sit empty for most subplots
        # TODO: Implement this with EdgeStack
        coltransform = mtransforms.blended_transform_factory(
            self.transAxes, self.figure.transFigure
        )
        rowtransform = mtransforms.blended_transform_factory(
            self.figure.transFigure, self.transAxes
        )
        self._left_label = self.text(
            0, 0.5, '', va='center', ha='right', transform=rowtransform
        )
        self._right_label = self.text(
            0, 0.5, '', va='center', ha='left', transform=rowtransform
        )
        self._bottom_label = self.text(
            0.5, 0, '', va='top', ha='center', transform=coltransform
        )
        self._top_label = self.text(
            0.5, 0, '', va='bottom', ha='center', transform=coltransform
        )

        # Axes inset title labels
        transform = self.transAxes
        self._upper_left_title = self.text(
            0, 0, '', va='top', ha='left', transform=transform,
        )
        self._upper_center_title = self.text(
            0, 0, '', va='top', ha='center', transform=transform,
        )
        self._upper_right_title = self.text(
            0, 0, '', va='top', ha='right', transform=transform,
        )
        self._lower_left_title = self.text(
            0, 0, '', va='bottom', ha='left', transform=transform,
        )
        self._lower_center_title = self.text(
            0, 0, '', va='bottom', ha='center', transform=transform,
        )
        self._lower_right_title = self.text(
            0, 0, '', va='bottom', ha='right', transform=transform,
        )

        # Abc label
        self._abc_label = self.text(0, 0, '', transform=transform)

        # Automatic axis sharing and formatting
        # TODO: Instead of format() call specific setters
        self._auto_share_setup()
        self.format(rc_mode=1)  # mode == 1 applies the rcShortParams

    def _auto_share_setup(self):
        """
        Automatically configure axis sharing based on the horizontal and
        vertical extent of subplots in the figure gridspec.
        """
        # Panel axes sharing, between main subplot and its panels
        # NOTE: _panel_share means "include this panel in the axis sharing group"
        # while _panel_sharex_group indicates the group itself and may include main axes
        def shared(paxs):
            return [pax for pax in paxs if not pax._panel_hidden and pax._panel_share]

        # Internal axis sharing, share stacks of panels and main axes with each other
        # NOTE: *This* block is why, even though share[xy] are figure-wide
        # settings, we still need the axes-specific _share[xy]_override attr
        if not self._panel_side:  # this is a main axes
            # Top and bottom
            bottom = self
            paxs = shared(self._bottom_panels)
            if paxs:
                bottom = paxs[-1]
                bottom._panel_sharex_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharex_group = True
                    iax._sharex_setup(bottom)  # parent is bottom-most
            paxs = shared(self._top_panels)
            for iax in paxs:
                iax._panel_sharex_group = True
                iax._sharex_setup(bottom)
            # Left and right
            # NOTE: Order of panel lists is always inside-to-outside
            left = self
            paxs = shared(self._left_panels)
            if paxs:
                left = paxs[-1]
                left._panel_sharey_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharey_group = True
                    iax._sharey_setup(left)  # parent is left-most
            paxs = shared(self._right_panels)
            for iax in paxs:
                iax._panel_sharey_group = True
                iax._sharey_setup(left)

        # External axes sharing, sometimes overrides panel axes sharing
        # NOTE: This can get very repetitive, but probably minimal impact?
        # Share x axes
        parent, *children = self._get_extent_axes('x')
        for child in children:
            child._sharex_setup(parent)
        # Share y axes
        parent, *children = self._get_extent_axes('y')
        for child in children:
            child._sharey_setup(parent)

    def _draw_auto_legends_colorbars(self):
        """
        Generate automatic legends and colorbars. Wrapper funcs
        let user add handles to location lists with successive calls to
        make successive calls to plotting commands.
        """
        for loc, (handles, kwargs) in self._auto_colorbar.items():
            self.colorbar(handles, **kwargs)
        for loc, (handles, kwargs) in self._auto_legend.items():
            self.legend(handles, **kwargs)
        self._auto_legend = {}
        self._auto_colorbar = {}

    def _get_extent_axes(self, x, panels=False):
        """
        Return the axes whose horizontal or vertical extent in the main
        gridspec matches the horizontal or vertical extend of this axes.
        The lefmost or bottommost axes are at the start of the list.
        """
        if not hasattr(self, 'get_subplotspec'):
            return [self]
        y = 'y' if x == 'x' else 'x'
        idx = 0 if x == 'x' else 1
        argfunc = np.argmax if x == 'x' else np.argmin
        irange = self._range_gridspec(x)
        if panels:
            axs = self.figure._iter_axes(hidden=False, children=False)
        else:
            axs = self.figure._axes_main
        axs = [ax for ax in axs if ax._range_gridspec(x) == irange]
        if not axs:
            return [self]
        else:
            pax = axs.pop(argfunc([ax._range_gridspec(y)[idx] for ax in axs]))
            return [pax, *axs]

    def _get_side_axes(self, side, panels=False):
        """
        Return the axes whose left, right, top, or bottom sides abutt
        against the same row or column as this axes.
        """
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid side {side!r}.')
        if not hasattr(self, 'get_subplotspec'):
            return [self]
        x = 'x' if side in ('left', 'right') else 'y'
        idx = 0 if side in ('left', 'top') else 1  # which side to test
        coord = self._range_gridspec(x)[idx]  # side for a particular axes
        if panels:
            axs = self.figure._iter_axes(hidden=False, children=False)
        else:
            axs = self.figure._axes_main
        axs = [ax for ax in axs if ax._range_gridspec(x)[idx] == coord]
        if not axs:
            return [self]
        else:
            return axs

    def _get_title(self, loc):
        """
        Get the title at the corresponding location.
        """
        if loc == 'abc':
            return self._abc_label
        else:
            return getattr(self, '_' + loc.replace(' ', '_') + '_title')

    def _hide_panel(self):
        """
        Hide axes contents but do *not* make the entire axes invisible. This
        is used to fill "panels" surreptitiously added to the gridspec
        for the purpose of drawing outer colorbars and legends.
        """
        # NOTE: Do not run self.clear in case we want to add a subplot title
        # above a colorbar on a top panel (see _reassign_title).
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        self._panel_hidden = True

    def _is_panel(self):
        """
        Return whether the current axes is a panel.
        """
        return bool(self._panel_parent)

    def _is_panel_parent_or_child(self, other):
        """
        Return whether the axes are related.
        """
        return self._panel_parent is other or other._panel_parent is self

    def _loc_translate(self, loc, mode=None, allow_manual=True):
        """
        Return the location string `loc` translated into a standardized form.
        """
        if mode == 'legend':
            valid = tuple(LOC_TRANSLATE.values())
        elif mode == 'panel':
            valid = ('left', 'right', 'top', 'bottom')
        elif mode == 'colorbar':
            valid = (
                'best', 'left', 'right', 'top', 'bottom',
                'upper left', 'upper right', 'lower left', 'lower right',
            )
        elif mode in ('abc', 'title'):
            valid = (
                'left', 'center', 'right',
                'upper left', 'upper center', 'upper right',
                'lower left', 'lower center', 'lower right',
            )
        else:
            raise ValueError(f'Invalid mode {mode!r}.')
        loc_translate = {
            key: value for key, value in LOC_TRANSLATE.items()
            if value in valid
        }
        if loc in (None, True):
            context = mode in ('abc', 'title')
            loc = rc.get(mode + '.loc', context=context)
            if loc is not None:
                loc = self._loc_translate(loc, mode)
        elif isinstance(loc, (str, Integral)):
            if loc in loc_translate.values():  # full name
                pass
            else:
                try:
                    loc = loc_translate[loc]
                except KeyError:
                    raise KeyError(
                        f'Invalid {mode} location {loc!r}. Options are: '
                        + ', '.join(map(repr, loc_translate)) + '.'
                    )
        elif (
            allow_manual
            and mode == 'legend'
            and np.iterable(loc)
            and len(loc) == 2
            and all(isinstance(l, Number) for l in loc)
        ):
            loc = np.array(loc)
        else:
            raise KeyError(f'Invalid {mode} location {loc!r}.')
        if mode == 'colorbar' and loc == 'best':  # white lie
            loc = 'lower right'
        return loc

    def _make_inset_locator(self, bounds, trans):
        """
        Return a locator that determines inset axes bounds.
        """
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
            return bb
        return inset_locator

    def _range_gridspec(self, x):
        """
        Return the column or row gridspec range for the axes.
        """
        if not hasattr(self, 'get_subplotspec'):
            raise RuntimeError('Axes is not a subplot.')
        ss = self.get_subplotspec()
        if hasattr(ss, 'get_active_rows_columns'):
            func = ss.get_active_rows_columns
        else:
            func = ss.get_rows_columns
        if x == 'x':
            _, _, _, _, col1, col2 = func()
            return col1, col2
        else:
            _, _, row1, row2, _, _ = func()
            return row1, row2

    def _range_tightbbox(self, x):
        """
        Return the tight bounding box span from the cached bounding box.
        `~proplot.axes.Axes.get_tightbbox` caches bounding boxes when
        `~Figure.get_tightbbox` is called.
        """
        # TODO: Better testing for axes visibility
        bbox = self._tightbbox
        if bbox is None:
            return np.nan, np.nan
        if x == 'x':
            return bbox.xmin, bbox.xmax
        else:
            return bbox.ymin, bbox.ymax

    def _reassign_subplot_label(self, side):
        """
        Re-assign the column and row labels to the relevant panel if
        present. This is called by `~proplot.figure.Figure._align_suplabel`.
        """
        # Place column and row labels on panels instead of axes -- works when
        # this is called on the main axes *or* on the relevant panel itself
        # TODO: Mixed figure panels with super labels? How does that work?
        if side == self._panel_side:
            ax = self._panel_parent
        else:
            ax = self
        paxs = getattr(ax, '_' + side + '_panels')
        if not paxs:
            return ax
        kw = {}
        pax = paxs[-1]  # outermost is always list in list
        obj = getattr(ax, '_' + side + '_label')
        for key in ('color', 'fontproperties'):  # TODO: add to this?
            kw[key] = getattr(obj, 'get_' + key)()
        pobj = getattr(pax, '_' + side + '_label')
        pobj.update(kw)
        text = obj.get_text()
        if text:
            obj.set_text('')
            pobj.set_text(text)
        return pax

    def _reassign_title(self):
        """
        Re-assign the title to the first upper panel if present. We cannot
        simply add the upper panel as a child axes, because then the title will
        be offset but still belong to main axes, which messes up the tight
        bounding box.
        """
        # Reassign title from main axes to top panel -- works when this is
        # called on the main axes *or* on the top panel itself. This is
        # critical for bounding box calcs; not always clear whether draw() and
        # get_tightbbox() are called on the main axes or panel first
        if self._panel_side == 'top' and self._panel_parent:
            ax = self._panel_parent
        else:
            ax = self
        taxs = ax._top_panels
        if not taxs or not ax._above_top_panels:
            tax = ax
        else:
            tax = taxs[-1]
            tax._title_pad = ax._title_pad
            for loc in ('abc', 'left', 'center', 'right'):
                kw = {}
                obj = ax._get_title(loc)
                if not obj.get_text():
                    continue
                tobj = tax._get_title(loc)
                for key in ('text', 'color', 'fontproperties'):  # add to this?
                    kw[key] = getattr(obj, 'get_' + key)()
                tobj.update(kw)
                obj.set_text('')

    def _sharex_setup(self, sharex):
        """
        Configure x-axis sharing for panels. Main axis sharing is done in
        `~CartesianAxes._sharex_setup`.
        """
        self._share_short_axis(sharex, 'left')  # x axis of left panels
        self._share_short_axis(sharex, 'right')
        self._share_long_axis(sharex, 'bottom')  # x axis of bottom panels
        self._share_long_axis(sharex, 'top')

    def _sharey_setup(self, sharey):
        """
        Configure y-axis sharing for panels. Main axis sharing is done in
        `~CartesianAxes._sharey_setup`.
        """
        self._share_short_axis(sharey, 'bottom')  # y axis of bottom panels
        self._share_short_axis(sharey, 'top')
        self._share_long_axis(sharey, 'left')  # y axis of left panels
        self._share_long_axis(sharey, 'right')

    def _share_short_axis(self, share, side):
        """
        Share the "short" axes of panels belonging to this subplot
        with panels belonging to an external subplot.
        """
        if share is None or self._panel_side:
            return  # if this is a panel
        axis = 'x' if side in ('left', 'right') else 'y'
        caxs = getattr(self, '_' + side + '_panels')
        paxs = getattr(share, '_' + side + '_panels')
        caxs = [pax for pax in caxs if not pax._panel_hidden]
        paxs = [pax for pax in paxs if not pax._panel_hidden]
        for cax, pax in zip(caxs, paxs):  # may be uneven
            getattr(cax, '_share' + axis + '_setup')(pax)

    def _share_long_axis(self, share, side):
        """
        Share the "long" axes of panels belonging to this subplot
        with panels belonging to an external subplot.
        """
        # NOTE: We do not check _panel_share because that only controls
        # sharing with main subplot, not other subplots
        if share is None or self._panel_side:
            return  # if this is a panel
        axis = 'x' if side in ('top', 'bottom') else 'y'
        paxs = getattr(self, '_' + side + '_panels')
        paxs = [pax for pax in paxs if not pax._panel_hidden]
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share)

    def _update_title_position(self, renderer):
        """
        Update the position of proplot inset titles and builtin matplotlib
        titles.
        """
        # Custom inset titles
        width, height = self.get_size_inches()
        for loc in (
            'abc',
            'upper left', 'upper right', 'upper center',
            'lower left', 'lower right', 'lower center',
        ):
            obj = self._get_title(loc)
            if loc == 'abc':
                loc = self._abc_loc
                if loc in ('left', 'right', 'center'):
                    continue
            if loc in ('upper center', 'lower center'):
                x = 0.5
            elif loc in ('upper left', 'lower left'):
                pad = rc['axes.titlepad'] / (72 * width)
                x = 1.5 * pad
            elif loc in ('upper right', 'lower right'):
                pad = rc['axes.titlepad'] / (72 * width)
                x = 1 - 1.5 * pad
            if loc in ('upper left', 'upper right', 'upper center'):
                pad = rc['axes.titlepad'] / (72 * height)
                y = 1 - 1.5 * pad
            elif loc in ('lower left', 'lower right', 'lower center'):
                pad = rc['axes.titlepad'] / (72 * height)
                y = 1.5 * pad
            obj.set_position((x, y))

        # Push title above tick marks, since builtin algorithm used to offset
        # the title seems to ignore them. This is known matplotlib problem but
        # especially annoying with top panels.
        # TODO: Make sure this is robust. Seems 'default' is returned usually
        # when tick sides is actually *both*.
        pad = self._title_pad
        pos = self.xaxis.get_ticks_position()
        fmt = self.xaxis.get_major_formatter()
        if self.xaxis.get_visible() and (
            pos == 'default'
            or (pos == 'top' and isinstance(fmt, mticker.NullFormatter))
            or (
                pos == 'unknown'
                and self._panel_side == 'top'
                and isinstance(fmt, mticker.NullFormatter)
            )
        ):
            pad += self.xaxis.get_tick_padding()
        pad_active = self._title_pad_active
        if pad_active is None or not np.isclose(pad_active, pad):
            # Avoid doing this on every draw in case it is expensive to change
            # the title Text transforms every time.
            self._title_pad_active = pad
            self._set_title_offset_trans(pad)

        # Adjust the title positions with builtin algorithm and match
        # the a-b-c text to the relevant position.
        super()._update_title_position(renderer)
        if self._abc_loc in ('left', 'center', 'right'):
            title = self._get_title(self._abc_loc)
            self._abc_label.set_position(title.get_position())
            self._abc_label.set_transform(
                self.transAxes + self.titleOffsetTrans
            )

    @staticmethod
    @warnings._rename_kwargs(mode='rc_mode')
    def _parse_format(rc_kw=None, rc_mode=None, **kwargs):
        """
        Separate `~proplot.config.rc` setting name value pairs from
        `~Axes.format` keyword arguments.
        """
        kw = {}
        rc_kw = rc_kw or {}
        rc_mode = _not_none(rc_mode, 2)
        for key, value in kwargs.items():
            key_fixed = rcsetup._rc_nodots.get(key, None)
            if key_fixed is None:
                kw[key] = value
            else:
                rc_kw[key_fixed] = value
        return rc_kw, rc_mode, kw

    def format(
        self, *, title=None, abovetop=None,
        figtitle=None, suptitle=None, rowlabels=None, collabels=None,
        leftlabels=None, rightlabels=None, toplabels=None, bottomlabels=None,
        llabels=None, rlabels=None, tlabels=None, blabels=None,
        ltitle=None, ctitle=None, rtitle=None,
        ultitle=None, uctitle=None, urtitle=None,
        lltitle=None, lctitle=None, lrtitle=None,
    ):
        """
        Modify the axes title(s), the a-b-c label, row and column labels, and
        the figure title. Called by the `~proplot.axes.CartesianAxes`,
        `~proplot.axes.PolarAxes`, and `~proplot.axes.GeoAxes` ``format``
        methods.

        Parameters
        ----------
        title : str, optional
            The axes title.
        abc : bool, optional
            Whether to apply "a-b-c" subplot labeling based on the subplot
            `~Axes.number`. If `~Axes.number` is greater than 26, the labels
            will loop around to a, ..., z, aa, ..., zz, aaa, ..., zzz, etc.
            Default is :rc:`abc`.
        abcstyle : str, optional
            String denoting the format of a-b-c labels containing the character
            ``a`` or ``A``. ``'a'`` is the default, but e.g. ``'a.'``,
            ``'a)'``, or ``'A'`` might also be desirable. Default is
            :rc:`abc.style`.
        abcloc, titleloc : str, optional
            Strings indicating the location for the a-b-c label and
            main title. The following locations keys are valid (defaults are
            :rc:`abc.loc` and :rc:`title.loc`):

            .. _title_table:

            ========================  ============================
            Location                  Valid keys
            ========================  ============================
            center above axes         ``'center'``, ``'c'``
            left above axes           ``'left'``, ``'l'``
            right above axes          ``'right'``, ``'r'``
            lower center inside axes  ``'lower center'``, ``'lc'``
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
        abovetop : bool, optional
            Whether to try to put the title and a-b-c label above the top panel
            (if it exists), or to always put them above the main subplot.
            Default is ``True``.
        ltitle, ctitle, rtitle, ultitle, uctitle, urtitle, lltitle, lctitle, \
lrtitle : str, optional
            Axes titles in specific positions (see `abcloc`). This lets you
            specify multiple title-like labels for a single subplot.
        leftlabels, rightlabels, toplabels, bottomlabels : list of str, \
optional
            Labels for the subplots lying along the left, right, top, and
            bottom edges of the figure. The length of each list must match
            the number of subplots along the corresponding edge.
        rowlabels, collabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`.
        llabels, rlabels, tlabels, blabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`, `rightlabels`,
            `bottomlabels`.
        figtitle, suptitle : str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.
            This is an improvement on matplotlib's "super" title, which just
            centers the text between figure edges.

        Note
        ----
        The `abc`, `abcstyle`, `abcloc`, and `titleloc` keyword arguments are
        actually :ref:`configuration settings <ug_config>` that are temporarily
        changed by the call to `~proplot.config.RcConfigurator.context`.  They
        are also documented here because it is very common to change them with
        `~Axes.format`.

        See also
        --------
        proplot.config.RcConfigurator.context
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        """
        # Misc axes settings
        # TODO: Add more settings to this?
        cycle = rc.get('axes.prop_cycle', context=True)
        if cycle is not None:
            self.set_prop_cycle(cycle)

        # Figure patch (for some reason needs to be re-asserted even if
        # declared before figure is drawn)
        kw = rc.fill({'facecolor': 'figure.facecolor'}, context=True)
        self.figure.patch.update(kw)
        if abovetop is not None:
            self._above_top_panels = abovetop
        pad = rc.get('axes.titlepad', context=True)
        if pad is not None:
            self._set_title_offset_trans(pad)
            self._title_pad = pad

        # Super title
        # NOTE: These are actually *figure-wide* settings, but that line
        # gets blurred where we have shared axes, spanning labels, and
        # whatnot. May result in redundant assignments if formatting more than
        # one axes, but operations are fast so some redundancy is nbd.
        # NOTE: Below kludge prevents changed *figure-wide* settings
        # from getting overwritten when user makes a new axes.
        fig = self.figure
        suptitle = _not_none(figtitle=figtitle, suptitle=suptitle)
        if len(fig._axes_main) > 1 and rc._context and rc._context[-1].mode == 1:
            kw = {}
        else:
            kw = rc.fill(
                {
                    'fontsize': 'suptitle.size',
                    'weight': 'suptitle.weight',
                    'color': 'suptitle.color',
                    'fontfamily': 'font.family'
                },
                context=True,
            )
        if suptitle or kw:
            fig._update_super_title(suptitle, **kw)

        # Labels
        rlabels = _not_none(rightlabels=rightlabels, rlabels=rlabels)
        blabels = _not_none(bottomlabels=bottomlabels, blabels=blabels)
        llabels = _not_none(
            rowlabels=rowlabels, leftlabels=leftlabels, llabels=llabels,
        )
        tlabels = _not_none(
            collabels=collabels, toplabels=toplabels, tlabels=tlabels,
        )
        for side, labels in zip(
            ('left', 'right', 'top', 'bottom'),
            (llabels, rlabels, tlabels, blabels)
        ):
            kw = rc.fill(
                {
                    'fontsize': side + 'label.size',
                    'weight': side + 'label.weight',
                    'color': side + 'label.color',
                    'fontfamily': 'font.family'
                },
                context=True,
            )
            if labels or kw:
                fig._update_subplot_labels(self, side, labels, **kw)

        # Helper function
        def sanitize_kw(kw, loc):
            kw = kw.copy()
            if loc in ('left', 'right', 'center'):
                kw.pop('border', None)
                kw.pop('borderwidth', None)
            return kw

        # A-b-c labels
        abc = False
        if not self._panel_side:
            # Properties
            # NOTE: Border props only apply for "inner" title locations so we
            # need to store on the axes whenever they are modified and always
            # re-apply the ones stored on the axes.
            kw = rc.fill(
                {
                    'fontsize': 'abc.size',
                    'weight': 'abc.weight',
                    'color': 'abc.color',
                    'fontfamily': 'font.family',
                },
                context=True
            )
            kwb = rc.fill(
                {
                    'border': 'abc.border',
                    'borderwidth': 'abc.borderwidth',
                },
                context=True,
            )
            self._abc_border_kwargs.update(kwb)
            kw.update(self._abc_border_kwargs)

            # Label format
            abcstyle = rc.get('abc.style', context=True)  # 1st run, or changed
            if abcstyle and self.number is not None:
                if not isinstance(abcstyle, str) or (
                    'a' not in abcstyle and 'A' not in abcstyle
                ):
                    raise ValueError(
                        f'Invalid abcstyle {abcstyle!r}. '
                        'Must include letter "a" or "A".'
                    )
                # Build abc labels as a...z...aa...zz...aaa...zzz
                # Permit abcstyles with arbitrary counts of A's
                nabc, iabc = divmod(self.number - 1, 26)
                text = (nabc + 1) * ABC_STRING[iabc]
                text = abcstyle.replace('a', text).replace('A', text.upper())
                self._abc_text = text

            # Apply text
            obj = self._abc_label
            abc = rc.get('abc', context=True)
            if abc is not None:
                obj.set_text(self._abc_text if bool(abc) else '')

            # Apply new settings
            loc = self._loc_translate(None, 'abc')
            loc_prev = self._abc_loc
            if loc is None:
                loc = loc_prev
            kw = sanitize_kw(kw, loc)
            if loc_prev is None or loc != loc_prev:
                obj_ref = self._get_title(loc)
                obj.set_ha(obj_ref.get_ha())
                obj.set_va(obj_ref.get_va())
                obj.set_transform(obj_ref.get_transform())
                obj.set_position(obj_ref.get_position())
            obj.update(kw)
            self._abc_loc = loc

        # Titles, with two workflows here:
        # 1. title='name' and titleloc='position'
        # 2. ltitle='name', rtitle='name', etc., arbitrarily many titles
        # NOTE: Matplotlib added axes.titlecolor in version 3.2 but we
        # still use custom title.size, title.weight, title.color
        # properties for retroactive support in older matplotlib versions.
        # First get params and update kwargs
        kw = rc.fill(
            {
                'fontsize': 'title.size',
                'weight': 'title.weight',
                'color': 'title.color',
                'fontfamily': 'font.family',
            },
            context=True
        )
        if 'color' in kw and kw['color'] == 'auto':
            del kw['color']  # WARNING: matplotlib permits invalid color here
        kwb = rc.fill(
            {
                'border': 'title.border',
                'borderwidth': 'title.borderwidth',
            },
            context=True,
        )
        self._title_border_kwargs.update(kwb)
        kw.update(self._title_border_kwargs)

        # Workflow 2, want this to come first so workflow 1 gets priority
        for iloc, ititle in zip(
            ('l', 'r', 'c', 'ul', 'uc', 'ur', 'll', 'lc', 'lr'),
            (
                ltitle, rtitle, ctitle,
                ultitle, uctitle, urtitle, lltitle, lctitle, lrtitle
            ),
        ):
            iloc = self._loc_translate(iloc, 'title')
            ikw = sanitize_kw(kw, iloc)
            iobj = self._get_title(iloc)
            iobj.update(ikw)
            if ititle is not None:
                iobj.set_text(ititle)

        # Workflow 1, make sure that if user calls ax.format(title='Title')
        # *then* ax.format(titleloc='left') it copies over the text.
        # Get current and previous location, prevent overwriting abc label
        loc = self._loc_translate(None, 'title')
        loc_prev = self._title_loc
        if loc is None:  # never None first run
            loc = loc_prev  # never None on subsequent runs

        # Reset old text
        if loc_prev is not None and loc != loc_prev:
            obj_prev = self._get_title(loc_prev)
            if title is None:
                title = obj_prev.get_text()
            obj_prev.set_text('')

        # Update new text
        kw = sanitize_kw(kw, loc)
        obj = self._get_title(loc)
        obj.update(kw)
        if title is not None:
            obj.set_text(title)
        self._title_loc = loc  # assigns default loc on first run

    def area(self, *args, **kwargs):
        """
        Alias for `~matplotlib.axes.Axes.fill_between`.
        """
        # NOTE: *Cannot* assign area = axes.Axes.fill_between because the
        # wrapper won't be applied and for some reason it messes up
        # automodsumm, which tries to put the matplotlib docstring on website
        return self.fill_between(*args, **kwargs)

    def areax(self, *args, **kwargs):
        """
        Alias for `~matplotlib.axes.Axes.fill_betweenx`.
        """
        return self.fill_betweenx(*args, **kwargs)

    def boxes(self, *args, **kwargs):
        """
        Alias for `~matplotlib.axes.Axes.boxplot`.
        """
        return self.boxplot(*args, **kwargs)

    def colorbar(
        self, *args, loc=None, pad=None,
        length=None, shrink=None, width=None, space=None, frame=None, frameon=None,
        alpha=None, linewidth=None, edgecolor=None, facecolor=None,
        **kwargs
    ):
        """
        Add an *inset* colorbar or *outer* colorbar along the outside edge of
        the axes. See `~proplot.axes.colorbar_wrapper` for details.

        Parameters
        ----------
        loc : str, optional
            The colorbar location. Default is :rc:`colorbar.loc`. The
            following location keys are valid:

            .. _colorbar_table:

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            outer left          ``'left'``, ``'l'``
            outer right         ``'right'``, ``'r'``
            outer bottom        ``'bottom'``, ``'b'``
            outer top           ``'top'``, ``'t'``
            default inset       ``'best'``, ``'inset'``, ``'i'``, ``0``
            upper right inset   ``'upper right'``, ``'ur'``, ``1``
            upper left inset    ``'upper left'``, ``'ul'``, ``2``
            lower left inset    ``'lower left'``, ``'ll'``, ``3``
            lower right inset   ``'lower right'``, ``'lr'``, ``4``
            "filled"            ``'fill'``
            ==================  =======================================

        pad : float or str, optional
            The space between the axes edge and the colorbar. For inset
            colorbars only. Units are interpreted by `~proplot.utils.units`.
            Default is :rc:`colorbar.insetpad`.
        length : float or str, optional
            The colorbar length. For outer colorbars, units are relative to the
            axes width or height. Default is :rc:`colorbar.length`. For inset
            colorbars, units are interpreted by `~proplot.utils.units`. Default
            is :rc:`colorbar.insetlength`.
        shrink : float, optional
            Alias for `length`. This is included to match the
            `matplotlib.figure.Figure.colorbar` keyword that has roughly the same
            meaning as `length`.
        width : float or str, optional
            The colorbar width. Units are interpreted by
            `~proplot.utils.units`.  For outer colorbars, default is
            :rc:`colorbar.width`. For inset colorbars, default is
            :rc:`colorbar.insetwidth`.
        space : float or str, optional
            For outer colorbars only. The space between the colorbar and the
            main axes. Units are interpreted by `~proplot.utils.units`.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, the default is :rc:`subplots.panelpad`.
        frame, frameon : bool, optional
            For inset colorbars, indicates whether to draw a "frame", just
            like `~matplotlib.axes.Axes.legend`. Default is
            :rc:`colorbar.frameon`.
        alpha, linewidth, edgecolor, facecolor : optional
            Transparency, edge width, edge color, and face color for the frame
            around the inset colorbar. Default is
            :rc:`colorbar.framealpha`, :rc:`axes.linewidth`,
            :rc:`axes.edgecolor`, and :rc:`axes.facecolor`,
            respectively.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~proplot.axes.colorbar_wrapper`.
        """
        # TODO: add option to pad inset away from axes edge!
        # TODO: get "best" colorbar location from legend algorithm.
        kwargs.update({'edgecolor': edgecolor, 'linewidth': linewidth})
        length = _not_none(length=length, shrink=shrink)
        if loc != 'fill':
            loc = self._loc_translate(loc, 'colorbar')

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            return ax.colorbar(loc='fill', *args, length=length, **kwargs)

        # Filled colorbar
        if loc == 'fill':
            # Hide content
            self._hide_panel()

            # Get subplotspec for colorbar axes
            side = self._panel_side
            length = _not_none(length, rc['colorbar.length'])
            subplotspec = self.get_subplotspec()
            if length <= 0 or length > 1:
                raise ValueError(
                    f'Panel colorbar length must satisfy 0 < length <= 1, '
                    f'got length={length!r}.'
                )
            if side in ('bottom', 'top'):
                gridspec = pgridspec._GridSpecFromSubplotSpec(
                    nrows=1, ncols=3, wspace=0,
                    subplot_spec=subplotspec,
                    width_ratios=((1 - length) / 2, length, (1 - length) / 2),
                )
                subplotspec = gridspec[1]
            else:
                gridspec = pgridspec._GridSpecFromSubplotSpec(
                    nrows=3, ncols=1, hspace=0,
                    subplot_spec=subplotspec,
                    height_ratios=((1 - length) / 2, length, (1 - length) / 2),
                )
                subplotspec = gridspec[1]

            # Draw colorbar axes
            with self.figure._context_authorize_add_subplot():
                ax = self.figure.add_subplot(subplotspec, projection='cartesian')  # noqa: E501
            self.add_child_axes(ax)

            # Location
            if side is None:  # manual
                orientation = kwargs.pop('orientation', None)
                if orientation == 'vertical':
                    side = 'left'
                else:
                    side = 'bottom'
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
            orientation_user = kwargs.get('orientation', None)
            if orientation_user and orientation_user != orientation:
                warnings._warn_proplot(
                    f'Overriding input orientation={orientation_user!r}.'
                )
            ticklocation = _not_none(
                ticklocation=kwargs.pop('ticklocation', None),
                tickloc=kwargs.pop('tickloc', None),
                default=ticklocation,
            )
            kwargs.update({
                'orientation': orientation,
                'ticklocation': ticklocation
            })

        # Inset colorbar
        else:
            # Default props
            cbwidth, cblength = width, length
            width, height = self.get_size_inches()
            extend = units(_not_none(
                kwargs.get('extendsize', None), rc['colorbar.insetextend']
            ))
            cbwidth = units(_not_none(
                cbwidth, rc['colorbar.insetwidth']
            )) / height
            cblength = units(_not_none(
                cblength, rc['colorbar.insetlength']
            )) / width
            pad = units(_not_none(pad, rc['colorbar.insetpad']))
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
                fbounds = (
                    1 - 2 * xpad - cblength,
                    1 - 2 * ypad - cbwidth - xspace
                )
            elif loc == 'upper left':
                bounds = (xpad, 1 - ypad - cbwidth)
                fbounds = (0, 1 - 2 * ypad - cbwidth - xspace)
            elif loc == 'lower left':
                bounds = (xpad, ypad + xspace)
                fbounds = (0, 0)
            else:
                bounds = (1 - xpad - cblength, ypad + xspace)
                fbounds = (1 - 2 * xpad - cblength, 0)
            bounds = (bounds[0], bounds[1], cblength, cbwidth)
            fbounds = (
                fbounds[0], fbounds[1],
                2 * xpad + cblength, 2 * ypad + cbwidth + xspace
            )

            # Make frame
            # NOTE: We do not allow shadow effects or fancy edges effect.
            # Also keep zorder same as with legend.
            frameon = _not_none(
                frame=frame, frameon=frameon, default=rc['colorbar.frameon'],
            )
            if frameon:
                xmin, ymin, width, height = fbounds
                patch = mpatches.Rectangle(
                    (xmin, ymin), width, height,
                    snap=True, zorder=4, transform=self.transAxes
                )
                alpha = _not_none(alpha, rc['colorbar.framealpha'])
                linewidth = _not_none(linewidth, rc['axes.linewidth'])
                edgecolor = _not_none(edgecolor, rc['axes.edgecolor'])
                facecolor = _not_none(facecolor, rc['axes.facecolor'])
                patch.update({
                    'alpha': alpha,
                    'linewidth': linewidth,
                    'edgecolor': edgecolor,
                    'facecolor': facecolor
                })
                self.add_artist(patch)

            # Make axes
            from .cartesian import CartesianAxes
            locator = self._make_inset_locator(bounds, self.transAxes)
            bbox = locator(None, None)
            ax = CartesianAxes(self.figure, bbox.bounds, zorder=5)
            ax.set_axes_locator(locator)
            self.add_child_axes(ax)

            # Default keyword args
            orient = kwargs.pop('orientation', None)
            if orient is not None and orient != 'horizontal':
                warnings._warn_proplot(
                    f'Orientation for inset colorbars must be horizontal, '
                    f'ignoring orient={orient!r}.'
                )
            ticklocation = kwargs.pop('tickloc', None)
            ticklocation = kwargs.pop('ticklocation', None) or ticklocation
            if ticklocation is not None and ticklocation != 'bottom':
                warnings._warn_proplot(
                    'Inset colorbars can only have ticks on the bottom.'
                )
            kwargs.update({
                'orientation': 'horizontal', 'ticklocation': 'bottom'
            })
            kwargs.setdefault('maxn', 5)
            kwargs.setdefault('extendsize', extend)

        # Generate colorbar
        return colorbar_wrapper(ax, *args, **kwargs)

    def legend(self, *args, loc=None, width=None, space=None, **kwargs):
        """
        Add an *inset* legend or *outer* legend along the edge of the axes.
        See `~proplot.axes.legend_wrapper` for details.

        Parameters
        ----------
        loc : int or str, optional
            The legend location. The following location keys are valid:

            .. _legend_table:

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            outer left          ``'left'``, ``'l'``
            outer right         ``'right'``, ``'r'``
            outer bottom        ``'bottom'``, ``'b'``
            outer top           ``'top'``, ``'t'``
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
            "filled"            ``'fill'``
            ==================  =======================================

        width : float or str, optional
            For outer legends only. The space allocated for the legend box.
            This does nothing if :rcraw:`tight` is ``True``. Units are
            interpreted by `~proplot.utils.units`.
        space : float or str, optional
            For outer legends only. The space between the axes and the legend
            box. Units are interpreted by `~proplot.utils.units`.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, the default is :rc:`subplots.panelpad`.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~proplot.axes.legend_wrapper`.
        """
        if loc != 'fill':
            loc = self._loc_translate(loc, 'legend')
        if isinstance(loc, np.ndarray):
            loc = loc.tolist()

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            return ax.legend(*args, loc='fill', **kwargs)

        # Fill
        if loc == 'fill':
            # Hide content
            self._hide_panel()

            # Try to make handles and stuff flush against the axes edge
            kwargs.setdefault('borderaxespad', 0)
            frameon = _not_none(
                kwargs.get('frame', None), kwargs.get('frameon', None),
                rc['legend.frameon']
            )
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
        # Perform extra post-processing steps
        self._reassign_title()
        super().draw(renderer, *args, **kwargs)

    def get_size_inches(self):
        # Return the width and height of the axes in inches.
        width, height = self.figure.get_size_inches()
        bbox = self.get_position()
        width = width * abs(bbox.width)
        height = height * abs(bbox.height)
        return np.array([width, height])

    def get_tightbbox(self, renderer, *args, **kwargs):
        # Perform extra post-processing steps and cache the bounding
        # box as an attribute.
        self._reassign_title()
        bbox = super().get_tightbbox(renderer, *args, **kwargs)
        self._tightbbox = bbox
        return bbox

    def heatmap(self, *args, aspect=None, **kwargs):
        """
        Pass all arguments to `~matplotlib.axes.Axes.pcolormesh` then apply
        settings that are suitable for heatmaps: square grid boxes by default,
        major ticks at the center of each grid box, no minor ticks, and no gridlines.

        Parameters
        ----------
        aspect : {'equal', 'auto'} or float, optional
            Controls the aspect ratio of the axes. The aspect is of particular
            relevance for heatmaps since it may distort the heatmap, i.e. a grid box
            will not be square. This parameter is a shortcut for explicitly calling
            `~matplotlib.axes.set_aspect`.

            The default is :rc:`image.heatmap`. The options are:

            - ``'equal'``: Ensures an aspect ratio of 1. Grid boxes will be square.
            - ``'auto'``: The axes is kept fixed and the aspect is adjusted so
              that the data fit in the axes. In general, this will result in non-square
              grid boxes.
        """
        obj = self.pcolormesh(*args, **kwargs)
        aspect = _not_none(aspect, rc['image.aspect'])
        xlocator = ylocator = None
        if hasattr(obj, '_coordinates'):
            coords = obj._coordinates
            coords = (coords[1:, ...] + coords[:-1, ...]) / 2
            coords = (coords[:, 1:, :] + coords[:, :-1, :]) / 2
            xlocator, ylocator = coords[0, :, 0], coords[:, 0, 1]
        self.format(
            aspect=aspect,
            xgrid=False, ygrid=False, xtickminor=False, ytickminor=False,
            xlocator=xlocator, ylocator=ylocator,
        )
        return obj

    @docstring.add_snippets
    def inset(self, *args, **kwargs):
        """
        %(axes.inset)s
        """
        return self.inset_axes(*args, **kwargs)

    @docstring.add_snippets
    def inset_axes(
        self, bounds, transform=None, zorder=4,
        zoom=None, zoom_kw=None,
        proj=None, proj_kw=None, projection=None, projection_kw=None, basemap=None,
        **kwargs
    ):
        """
        %(axes.inset)s
        """
        # Carbon copy with my custom axes
        if not transform:
            transform = self.transAxes
        else:
            transform = _get_transform(self, transform)
        label = kwargs.pop('label', 'inset_axes')
        proj = _not_none(proj=proj, projection=projection)
        proj_kw = _not_none(proj_kw=proj_kw, projection_kw=projection_kw, default={})

        # Inherit from current axes
        if proj is None:
            proj = self.name
            if basemap is not None:
                proj_kw['basemap'] = basemap
            if proj_kw:
                warnings._warn_proplot(
                    'Inheriting projection from the main axes. '
                    f'Ignoring proj_kw keyword args: {proj_kw}'
                )
            if proj in ('cartopy', 'basemap'):
                map_projection = copy.copy(self.projection)
                kwargs.setdefault('map_projection', map_projection)

        # Create new projection
        elif proj == 'cartesian':
            pass
        elif proj == 'polar':
            proj = 'polar2'  # custom proplot name
        else:
            proj_kw.setdefault('basemap', basemap)
            map_projection = constructor.Proj(proj, **proj_kw)
            kwargs.setdefault('map_projection', map_projection)
            proj = 'basemap' if proj_kw['basemap'] else 'cartopy'

        # This puts the rectangle into figure-relative coordinates.
        locator = self._make_inset_locator(bounds, transform)
        cls = mprojections.get_projection_class(proj)
        bb = locator(None, None)
        ax = cls(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)

        # The following locator lets the axes move if we used data coordinates,
        # is called by ax.apply_aspect()
        zoom = _not_none(zoom, self.name == ax.name)  # only zoom when same projection
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        ax._inset_zoom = zoom
        ax._inset_parent = self

        # Zoom indicator (NOTE: Requires version >=3.0)
        if zoom:
            zoom_kw = zoom_kw or {}
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def indicate_inset_zoom(
        self, alpha=None,
        lw=None, linewidth=None, zorder=3.5,
        ls=None, linestyle=None,
        color=None, edgecolor=None, **kwargs
    ):
        """
        Draw lines indicating the zoom range of the inset axes. This is similar
        to the builtin `~matplotlib.axes.Axes.indicate_inset_zoom` except
        lines are *refreshed* at draw-time. This is also called automatically
        when ``zoom=True`` is passed to `~Axes.inset_axes`. Note this method
        must be called from the *inset* axes and not the parent axes.

        Parameters
        ----------
        alpha : float, optional
            The transparency of the zoom box fill.
        lw, linewidth : float, optional
            The width of the zoom lines and box outline in points.
        ls, linestyle : linestyle-spec, optional
            The line style for the zoom lines and box outline.
        color, edgecolor : color-spec, optional
            The color of the zoom lines and box outline.
        capstyle : {'butt', 'round', 'projecting'}
            The cap style for the zoom lines and box outline.
        zorder : float, optional
            The `zorder \
<https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html>`__
            of the axes, should be greater than the zorder of
            elements in the parent axes. Default is ``3.5``.

        Other parameters
        ----------------
        **kwargs
            Passed to `~matplotlib.axes.Axes.indicate_inset`.
        """
        # Should be called from the inset axes
        parent = self._inset_parent
        alpha = alpha or 1.0
        kwargs.setdefault('capstyle', 'round')  # match zoom capstyle
        linestyle = _not_none(ls=ls, linestyle=linestyle, default='-')
        linewidth = _not_none(
            lw=lw, linewidth=linewidth, default=rc['axes.linewidth'],
        )
        edgecolor = _not_none(
            color=color, edgecolor=edgecolor, default=rc['axes.edgecolor'],
        )
        if not parent:
            raise ValueError(f'{self} is not an inset axes.')
        xlim, ylim = self.get_xlim(), self.get_ylim()
        rect = (xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0])

        # Call indicate_inset
        props = {
            'linestyle': linestyle, 'linewidth': linewidth,
            'edgecolor': edgecolor, 'alpha': alpha,
        }
        rectpatch, connects = parent.indicate_inset(
            rect, self, zorder=zorder, **props, **kwargs
        )

        # Update zoom or adopt properties from old one
        if self._inset_zoom_data is None:
            for line in connects:
                line.update(props)
        else:
            rectpatch_prev, connects_prev = self._inset_zoom_data
            rectpatch.update_from(rectpatch_prev)
            rectpatch_prev.set_visible(False)
            for line, line_prev in zip(connects, connects_prev):
                visible = line.get_visible()
                line.update_from(line_prev)
                line.set_visible(visible)
                line_prev.set_visible(False)
        self._inset_zoom_data = (rectpatch, connects)
        return rectpatch, connects

    @docstring.add_snippets
    def panel(self, side, **kwargs):
        """
        %(axes.panel)s
        """
        return self.panel_axes(side, **kwargs)

    @docstring.add_snippets
    def panel_axes(self, side, **kwargs):
        """
        %(axes.panel)s
        """
        side = self._loc_translate(side, 'panel')
        return self.figure._add_axes_panel(self, side, **kwargs)

    @_parametric_wrapper
    @_standardize_1d
    @_cmap_changer
    def parametric(
        self, *args, values=None,
        cmap=None, norm=None, interp=0,
        scalex=True, scaley=True,
        **kwargs
    ):
        """
        Draw a line whose color changes as a function of the parametric
        coordinate ``values`` using the input colormap ``cmap``.
        Invoked when you pass the `cmap` keyword argument to
        `~matplotlib.axes.Axes.plot`.

        Parameters
        ----------
        *args : (y,), (x, y), or (x, y, values)
            The coordinates. If `x` is not provided, it is inferred from `y`.
        values : list of float
            The parametric values used to map points on the line to colors
            in the colormap. This can also be passed as a third positional argument.
        cmap : colormap spec, optional
            The colormap specifier, passed to `~proplot.constructor.Colormap`.
        cmap_kw : dict, optional
            Keyword arguments passed to `~proplot.constructor.Colormap`.
        norm : normalizer spec, optional
            The normalizer, passed to `~proplot.constructor.Norm`.
        norm_kw : dict, optional
            Keyword arguments passed to `~proplot.constructor.Norm`.
        interp : int, optional
            If greater than ``0``, we interpolate to additional points
            between the `values` coordinates. The number corresponds to the
            number of additional color levels between the line joints
            and the halfway points between line joints.
        scalex, scaley : bool, optional
            Whether the view limits are adapted to the data limits. The values are
            passed on to `~matplotlib.axes.Axes.autoscale_view`.

        Other parameters
        ----------------
        **kwargs
            Valid `~matplotlib.collections.LineCollection` properties.

        Returns
        -------
        `~matplotlib.collections.LineCollection`
            The parametric line. See `this matplotlib example \
<https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line>`__.
        """
        # Get x/y coordinates and values for points to the 'left' and 'right'
        # of each joint
        x, y = args  # standardized by parametric wrapper
        interp  # avoid U100 unused argument error (arg is handled by wrapper)
        coords = []
        levels = edges(values)
        for j in range(y.shape[0]):
            if j == 0:
                xleft, yleft = [], []
            else:
                xleft = [0.5 * (x[j - 1] + x[j]), x[j]]
                yleft = [0.5 * (y[j - 1] + y[j]), y[j]]
            if j + 1 == y.shape[0]:
                xright, yright = [], []
            else:
                xleft = xleft[:-1]  # prevent repetition when joined with right
                yleft = yleft[:-1]
                xright = [x[j], 0.5 * (x[j + 1] + x[j])]
                yright = [y[j], 0.5 * (y[j + 1] + y[j])]
            pleft = np.stack((xleft, yleft), axis=1)
            pright = np.stack((xright, yright), axis=1)
            coords.append(np.concatenate((pleft, pright), axis=0))
        coords = np.array(coords)

        # Create LineCollection and update with values
        # NOTE: Default capstyle is butt but this may look weird with vector graphics
        hs = mcollections.LineCollection(
            coords, cmap=cmap, norm=norm,
            linestyles='-', capstyle='butt', joinstyle='miter',
        )
        values = np.asarray(values)
        hs.set_array(values)
        hs.update({
            key: value for key, value in kwargs.items()
            if key not in ('color',)
        })

        # Add collection with some custom attributes
        # NOTE: Modern API uses self._request_autoscale_view but this is
        # backwards compatible to earliest matplotlib versions.
        self.add_collection(hs)
        self.autoscale_view(scalex=scalex, scaley=scaley)
        hs.values = values
        hs.levels = levels  # needed for other functions
        return hs

    def violins(self, *args, **kwargs):
        """
        Alias for `~matplotlib.axes.Axes.violinplot`.
        """
        return self.violinplot(*args, **kwargs)

    def _iter_axes(self, panels=None, hidden=False, children=False):
        """
        Return a list of axes and child panel axes.

        Parameters
        ----------
        panels : tuple, optional
            Tuple of panels to select.
        hidden : bool, optional
            Whether to include "hidden" panels.
        children : bool, optional
            Whether to include children.
        """
        panels = _not_none(panels, ('left', 'right', 'bottom', 'top'))
        if not set(panels) <= {'left', 'right', 'bottom', 'top'}:
            raise ValueError(f'Invalid sides {panels!r}.')
        for iax in (
            self,
            *(
                jax for side in panels
                for jax in getattr(self, '_' + side + '_panels')
            )
        ):
            if not hidden and iax._panel_hidden:
                continue  # ignore hidden panel and its colorbar/legend child
            for jax in ((iax, *iax.child_axes) if children else (iax,)):
                if not jax.get_visible():
                    continue  # safety first
                yield jax

    @property
    def number(self):
        """
        The axes number. This controls the order of a-b-c labels and the
        order of appearence in the `~proplot.ui.SubplotsContainer` returned by
        `~proplot.ui.subplots`.
        """
        return self._number

    @number.setter
    def number(self, num):
        if num is not None and (not isinstance(num, Integral) or num < 1):
            raise ValueError(f'Invalid number {num!r}. Must be integer >=1.')
        self._number = num

    # For consistency with _left_title, _upper_left_title, etc.
    _center_title = property(lambda self: self.title)

    # Wrapped by special functions
    # Also support redirecting to Basemap methods
    text = _text_wrapper(
        maxes.Axes.text
    )
    plot = _plot_wrapper(_standardize_1d(_indicate_error(_cycle_changer(
        maxes.Axes.plot
    ))))
    scatter = _scatter_wrapper(_standardize_1d(_indicate_error(_cycle_changer(
        maxes.Axes.scatter
    ))))
    bar = _bar_wrapper(_standardize_1d(_indicate_error(_cycle_changer(
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
    violinplot = _violinplot_wrapper(_standardize_1d(_indicate_error(
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
    step = _standardize_1d(_cycle_changer(
        maxes.Axes.step
    ))

    # Wrapped by standardizer
    stem = _standardize_1d(_stem_wrapper(
        maxes.Axes.stem
    ))
    hlines = _standardize_1d(_hlines_wrapper(
        maxes.Axes.hlines
    ))
    vlines = _standardize_1d(_vlines_wrapper(
        maxes.Axes.vlines
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
    pcolorfast = _standardize_2d(_cmap_changer(
        maxes.Axes.pcolorfast  # WARNING: not available in cartopy and basemap
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
