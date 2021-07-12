#!/usr/bin/env python3
"""
The base axes class used for all ProPlot figures.
"""
import copy
import re
from numbers import Integral, Number

import matplotlib.axes as maxes
import matplotlib.collections as mcollections
import matplotlib.legend as mlegend
import matplotlib.patches as mpatches
import matplotlib.projections as mprojections
import matplotlib.transforms as mtransforms
import numpy as np

from .. import constructor
from .. import gridspec as pgridspec
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import _not_none, docstring, rcsetup, warnings
from ..utils import edges, units
from . import plot as wrap

__all__ = ['Axes']

ABC_STRING = 'abcdefghijklmnopqrstuvwxyz'
KEYS_INNER = (
    'border', 'borderwidth', 'bbox', 'bboxpad', 'bboxcolor', 'bboxstyle', 'bboxalpha',
)
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
    Keyword arguments used to update the background patch. This can
    be used e.g. to apply background hatching with ``patch_kw={'hatch': 'xxx'}``.
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
transform : {'data', 'axes', 'figure'} or `~matplotlib.transforms.Transform`, optional
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
    The `zorder <https://matplotlib.org/stable/gallery/misc/zorder_demo.html>`__
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
    def __init__(self, *args, number=None, main=False, _subplotspec=None, **kwargs):
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
        # TODO: Why are some of these not set in __init__?
        if main:
            self.figure._subplots_main.append(self)
        self.number = number  # for a-b-c numbering
        self._auto_format = None  # manipulated by wrapper functions
        self._abc_loc = None
        self._abc_text = None
        self._abc_border_kwargs = {}
        self._abc_pad = rc['abc.titlepad']
        self._title_loc = None
        self._title_border_kwargs = {}  # title border properties
        self._title_above = rc['title.above']
        self._title_pad = rc['title.pad']
        self._title_pad_current = None
        self._tight_bbox = None  # bounding boxes are saved
        self._panel_hidden = False  # True when "filled" with cbar/legend
        self._panel_parent = None
        self._panel_share = False
        self._panel_sharex_group = False
        self._panel_sharey_group = False
        self._panel_side = None
        self._inset_parent = None
        self._inset_zoom = False
        self._inset_zoom_data = None

        # Axes colorbars and legends
        self._colorbar_dict = {}
        self._legend_dict = {}

        # Axes panels
        d = self._panel_dict = {}
        d['left'] = []  # NOTE: panels will be sorted inside-to-outside
        d['right'] = []
        d['bottom'] = []
        d['top'] = []

        # Axes titles
        # Record the original positions to account for offsetting
        d = self._title_dict = {}
        ta = self.transAxes
        d['abc'] = self.text(0, 0, '', transform=ta)
        d['left'] = self._left_title  # WARNING: track in case mpl changes this
        d['center'] = self.title
        d['right'] = self._right_title
        d['upper left'] = self.text(0, 0, '', va='top', ha='left', transform=ta)
        d['upper center'] = self.text(0, 0, '', va='top', ha='center', transform=ta)
        d['upper right'] = self.text(0, 0, '', va='top', ha='right', transform=ta)
        d['lower left'] = self.text(0, 0, '', va='bottom', ha='left', transform=ta)
        d['lower center'] = self.text(0, 0, '', va='bottom', ha='center', transform=ta)
        d['lower right'] = self.text(0, 0, '', va='bottom', ha='right', transform=ta)

        # Axes row and column labels
        # NOTE: Most of these sit empty for most subplots
        # TODO: Implement this with EdgeStack, avoid creating silly empty objects
        d = self._label_dict = {}
        tf = self.figure.transFigure
        tc = mtransforms.blended_transform_factory(ta, tf)
        tr = mtransforms.blended_transform_factory(tf, ta)
        d['left'] = self.text(0, 0.5, '', va='center', ha='right', transform=tr)
        d['right'] = self.text(0, 0.5, '', va='center', ha='left', transform=tr)
        d['bottom'] = self.text(0.5, 0, '', va='top', ha='center', transform=tc)
        d['top'] = self.text(0.5, 0, '', va='bottom', ha='center', transform=tc)
        d = self._label_pad = {}
        d['left'] = rc['leftlabel.pad']
        d['right'] = rc['rightlabel.pad']
        d['bottom'] = rc['bottomlabel.pad']
        d['top'] = rc['toplabel.pad']

        # Subplot spec
        # WARNING: For mpl>=3.4.0 subplotspec assigned *after* initialization using
        # set_subplotspec. Tried to defer to setter but really messes up both format()
        # and _auto_share_setup(). Instead use workaround: Have Figure.add_subplot pass
        # subplotspec as a hidden keyword arg. Non-subplots don't need this arg.
        # See https://github.com/matplotlib/matplotlib/pull/18564
        if _subplotspec is not None:
            self.set_subplotspec(_subplotspec)

        # Default sharing and formatting
        # TODO: Apply specific setters instead of format()
        self._auto_share_setup()
        self.format(rc_mode=1)  # rc_mode == 1 applies the custom proplot params

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
            paxs = shared(self._panel_dict['bottom'])
            if paxs:
                bottom = paxs[-1]
                bottom._panel_sharex_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharex_group = True
                    iax._sharex_setup(bottom)  # parent is bottom-most
            paxs = shared(self._panel_dict['top'])
            for iax in paxs:
                iax._panel_sharex_group = True
                iax._sharex_setup(bottom)
            # Left and right
            # NOTE: Order of panel lists is always inside-to-outside
            left = self
            paxs = shared(self._panel_dict['left'])
            if paxs:
                left = paxs[-1]
                left._panel_sharey_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharey_group = True
                    iax._sharey_setup(left)  # parent is left-most
            paxs = shared(self._panel_dict['right'])
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
            axs = self.figure._subplots_main
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
            axs = self.figure._subplots_main
        axs = [ax for ax in axs if ax._range_gridspec(x)[idx] == coord]
        if not axs:
            return [self]
        else:
            return axs

    def _hide_panel(self):
        """
        Hide axes contents but do *not* make the entire axes invisible. This is used to
        fill "panels" surreptitiously added to the gridspec for the purpose of drawing
        outer colorbars and legends.
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

    def _is_panel_group_member(self, other):
        """
        Return whether the axes are related.
        """
        return (
            self._panel_parent is other  # child
            or other._panel_parent is self   # parent
            or other._panel_parent is self._panel_parent  # sibling
        )

    def _loc_translate(self, loc, mode=None):
        """
        Return the location string `loc` translated into a standardized form.
        """
        if mode == 'legend':
            options = tuple(LOC_TRANSLATE.values())
        elif mode == 'panel':
            options = ('left', 'right', 'top', 'bottom')
        elif mode == 'colorbar':
            options = (
                'best', 'left', 'right', 'top', 'bottom',
                'upper left', 'upper right', 'lower left', 'lower right',
            )
        elif mode in ('abc', 'title'):
            options = (
                'left', 'center', 'right',
                'upper left', 'upper center', 'upper right',
                'lower left', 'lower center', 'lower right',
            )
        else:
            raise ValueError(f'Invalid mode {mode!r}.')
        loc_translate = {
            key: value
            for short, long in LOC_TRANSLATE.items()
            for key, value in ((long, long), (short, long))
            if long in options
        }
        if loc in (None, True):
            context = mode in ('abc', 'title')
            loc = rc.get(mode + '.loc', context=context)
            if loc is not None:
                loc = self._loc_translate(loc, mode)
        elif isinstance(loc, (str, Integral)):
            try:
                loc = loc_translate[loc]
            except KeyError:
                raise KeyError(
                    f'Invalid {mode} location {loc!r}. Options are: '
                    + ', '.join(map(repr, loc_translate)) + '.'
                )
        elif (
            mode == 'legend'
            and np.iterable(loc)
            and len(loc) == 2
            and all(isinstance(l, Number) for l in loc)
        ):
            loc = tuple(loc)
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
        bbox = self._tight_bbox
        if bbox is None:
            return np.nan, np.nan
        if x == 'x':
            return bbox.xmin, bbox.xmax
        else:
            return bbox.ymin, bbox.ymax

    def _reassign_label(self, side):
        """
        Reassign the column and row labels to the relevant panel if present.
        This is called by `~proplot.figure.Figure._align_subplot_figure_labels`.
        """
        # NOTE: Since panel axes are "children" main axes is always drawn first.
        paxs = self._panel_dict[side]
        if not paxs:
            return self
        kw = {}
        pax = paxs[-1]  # outermost
        cobj = self._label_dict[side]
        pobj = pax._label_dict[side]
        for key in ('text', 'color', 'fontproperties'):
            kw[key] = getattr(cobj, 'get_' + key)()
        pobj.update(kw)
        cobj.set_text('')
        return pax

    def _reassign_title(self):
        """
        Re-assign the title to the first upper panel if present. We cannot
        simply add the upper panel as a child axes, because then the title will
        be offset but still belong to main axes, which messes up the tight
        bounding box.
        """
        # NOTE: Since panel axes are "children" main axes is always drawn first.
        taxs = self._panel_dict['top']
        if not taxs or not self._title_above:
            return
        tax = taxs[-1]  # outermost
        tax._title_pad = self._title_pad
        for loc in ('abc', 'left', 'center', 'right'):
            kw = {}
            cobj = self._title_dict[loc]
            if loc == 'abc':
                loc = tax._abc_loc = self._abc_loc
            if loc not in ('left', 'center', 'right'):
                continue
            tobj = tax._title_dict[loc]
            text = cobj.get_text()
            if not text:
                continue
            for key in ('color', 'fontproperties'):
                kw[key] = getattr(cobj, 'get_' + key)()
            tobj.update(kw)
            tobj.set_text(text)
            cobj.set_text('')

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
        caxs = self._panel_dict[side]
        paxs = share._panel_dict[side]
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
        paxs = self._panel_dict[side]
        paxs = [pax for pax in paxs if not pax._panel_hidden]
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share)

    def _update_abc(self):
        """
        Whether to update the label.
        """
        abc = False
        if self._panel_side:
            return

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
                'bbox': 'abc.bbox',
                'bboxpad': 'abc.bboxpad',
                'bboxcolor': 'abc.bboxcolor',
                'bboxstyle': 'abc.bboxstyle',
                'bboxalpha': 'abc.bboxalpha',
            },
            context=True,
        )
        self._abc_border_kwargs.update(kwb)
        kw.update(self._abc_border_kwargs)

        # A-b-c labels. Build as a...z...aa...zz...aaa...zzz
        style = rc.get('abc.style', context=True)  # 1st run, or changed
        if style and self.number is not None:
            if not isinstance(style, str) or 'a' not in style and 'A' not in style:
                raise ValueError(
                    f'Invalid abcstyle {style!r}. Must include letter "a" or "A".'
                )
            nabc, iabc = divmod(self.number - 1, 26)
            old = re.search('[aA]', style).group()  # return the *first* 'a'
            new = (nabc + 1) * ABC_STRING[iabc]
            new = new.upper() if old == 'A' else new
            self._abc_text = style.replace(old, new, 1)

        # Apply a-b-c text
        abc = rc.get('abc', context=True)
        aobj = self._title_dict['abc']
        if abc is not None:
            aobj.set_text(self._abc_text if bool(abc) else '')

        # Apply a-b-c settings
        loc = self._loc_translate(None, 'abc')
        loc_prev = self._abc_loc
        if loc is None:
            loc = loc_prev
        if loc in ('left', 'right', 'center'):
            for key in KEYS_INNER:
                kw.pop(key, None)
        aobj.update(kw)
        self._abc_loc = loc

    def _update_super(self, suptitle, **kwargs):
        """
        Update super title and row and column labels.
        """
        # NOTE: These are actually *figure-wide* settings, but that line gets
        # blurred where we have shared axes, spanning labels, and whatnot. May result
        # in redundant assignments if formatting more than one axes, but operations
        # are fast so some redundancy is nbd.
        # NOTE: Below kludge prevents changed *figure-wide* settings from getting
        # overwritten when user makes a new panels or insets. Funky limnitation but
        # kind of makes sense if these are inaccessible from panels.
        fig = self.figure
        ignore = self not in fig._subplots_main
        kw = {} if ignore else rc.fill(
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
        for side, labels in kwargs.items():
            kw = {} if ignore else rc.fill(
                {
                    'fontsize': side + 'label.size',
                    'weight': side + 'label.weight',
                    'color': side + 'label.color',
                    'fontfamily': 'font.family'
                },
                context=True,
            )
            if labels or kw:
                fig._update_super_labels(self, side, labels, **kw)

    def _update_title_all(self, title=None, **kwargs):
        """
        Update the titles.
        """
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
                'bbox': 'title.bbox',
                'bboxpad': 'title.bboxpad',
                'bboxcolor': 'title.bboxcolor',
                'bboxstyle': 'title.bboxstyle',
                'bboxalpha': 'title.bboxalpha',
            },
            context=True,
        )
        self._title_border_kwargs.update(kwb)
        kw.update(self._title_border_kwargs)

        # Workflow 2, want this to come first so workflow 1 gets priority
        for iloc, ititle in kwargs.items():
            ikw = kw.copy()
            iloc = self._loc_translate(iloc, 'title')
            if iloc in ('left', 'center', 'right'):
                for key in KEYS_INNER:
                    ikw.pop(key, None)
            iobj = self._title_dict[iloc]
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

        # Remove previous text
        if loc_prev is not None and loc != loc_prev:
            tobj_prev = self._title_dict[loc_prev]
            if title is None:
                title = tobj_prev.get_text()
            tobj_prev.set_text('')

        # Add new text and settings
        kw = kw.copy()
        if loc in ('left', 'center', 'right'):
            for key in KEYS_INNER:
                kw.pop(key, None)
        tobj = self._title_dict[loc]
        tobj.update(kw)
        if title is not None:
            tobj.set_text(title)
        self._title_loc = loc  # assigns default loc on first run

    def _update_title_position(self, renderer):
        """
        Update the position of proplot inset titles and builtin matplotlib
        titles. This is called by matplotlib at drawtime.
        """
        # Update title positions
        # NOTE: Critical to do this every time in case padding changes or
        # we added or removed an a-b-c label in the same position as a title
        width, height = self.get_size_inches()
        x_pad = self._title_pad / (72 * width)
        y_pad = self._title_pad / (72 * height)
        for loc, obj in self._title_dict.items():
            x, y = (0, 1)
            if loc == 'abc':  # redirect
                loc = self._abc_loc
            if loc == 'left':
                x = 0
            elif loc == 'center':
                x = 0.5
            elif loc == 'right':
                x = 1
            if loc in ('upper center', 'lower center'):
                x = 0.5
            elif loc in ('upper left', 'lower left'):
                x = x_pad
            elif loc in ('upper right', 'lower right'):
                x = 1 - x_pad
            if loc in ('upper left', 'upper right', 'upper center'):
                y = 1 - y_pad
            elif loc in ('lower left', 'lower right', 'lower center'):
                y = y_pad
            obj.set_position((x, y))

        # Push title above tick marks, since builtin algorithm seems to ignore them.
        # This is known matplotlib problem but especially annoying with top panels.
        # NOTE: See axis.get_ticks_position for inspiration
        pad = self._title_pad
        if self.xaxis.get_visible() and any(
            tick.tick2line.get_visible() and not tick.label2.get_visible()
            for tick in self.xaxis.majorTicks
        ):
            pad += self.xaxis.get_tick_padding()

        # Avoid applying padding on every draw in case it is expensive to change
        # the title Text transforms every time.
        pad_current = self._title_pad_current
        if pad_current is None or not np.isclose(pad, pad_current):
            self._title_pad_current = pad
            self._set_title_offset_trans(pad)

        # Adjust the above-axes positions with builtin algorithm
        # WARNING: Make sure the name of this private function doesn't change
        super()._update_title_position(renderer)

        # Sync the title positiona with the a-b-c label position
        aobj = self._title_dict['abc']
        tobj = self._title_dict[self._abc_loc]
        aobj.set_ha(tobj.get_ha())
        aobj.set_va(tobj.get_va())
        aobj.set_position(tobj.get_position())
        aobj.set_transform(tobj.get_transform())

        # Offset title away from a-b-c label
        # NOTE: Title texts all use axes transform in x-direction
        # TODO: Make empirical padding of '0.4' em tunable?
        if not tobj.get_text() or not aobj.get_text():
            return
        awidth, twidth = (
            obj.get_window_extent(renderer).transformed(self.transAxes.inverted()).width
            for obj in (aobj, tobj)
        )
        apad, tpad = (
            (self._abc_pad / 72) / self.get_size_inches()[0]
            for obj in (aobj, tobj)
        )
        ha = aobj.get_ha()
        aoffset = toffset = 0
        if ha == 'left':
            toffset = awidth + apad
        elif ha == 'right':
            aoffset = -(twidth + tpad)
        else:  # guaranteed center, there are others
            toffset = 0.5 * (awidth + apad)
            aoffset = -0.5 * (twidth + tpad)
        aobj.set_x(aobj.get_position()[0] + aoffset)
        tobj.set_x(tobj.get_position()[0] + toffset)

    @staticmethod
    @warnings._rename_kwargs('0.6', mode='rc_mode')
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

    @docstring.add_snippets
    def format(
        self, *, title=None,
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

        ltitle, ctitle, rtitle, ultitle, uctitle, urtitle, lltitle, lctitle, lrtitle \
: str, optional
            Axes titles in specific positions. Works as an alternative to
            ``ax.format(title='title', titleloc='loc')`` and lets you specify
            multiple title-like labels in a single subplot.
        abcborder, titleborder : bool, optional
            Whether to draw a white border around titles and a-b-c labels
            positioned inside the axes. This can help them stand out on top
            of artists plotted inside the axes. Defaults are
            :rc:`abc.border` and :rc:`title.border`
        abcbbox, titlebbox : bool, optional
            Whether to draw a white bbox around titles and a-b-c labels
            positioned inside the axes. This can help them stand out on top
            of artists plotted inside the axes. Defaults are
            :rc:`abc.bbox` and :rc:`title.bbox`
        titlepad : float, optional
            The padding for the inner and outer titles and a-b-c labels in
            arbitrary units (default is points). Default is :rc:`title.pad`.
        titleabove : bool, optional
            Whether to try to put outer titles and a-b-c labels above panels,
            colorbars, or legends that are above the axes. Default is :rc:`title.above`.
        leftlabels, toplabels, rightlabels, bottomlabels : list of str, optional
            Labels for the subplots lying along the left, top, right, and
            bottom edges of the figure. The length of each list must match
            the number of subplots along the corresponding edge.
        rowlabels, collabels, llabels, tlabels, rlabels, blabels : list of str, optional
            Aliases for `leftlabels`, `toplabels`, `leftlabels`, `toplabels`,
            `rightlabels`, and `bottomlabels`.
        leftlabelpad, toplabelpad, rightlabelpad, bottomlabelpad : float, optional
            The padding between the labels and the axes content in arbitrary units
            (default is points). Defaults are :rcraw:`leftlabel.pad`,
            :rcraw:`toplabel.pad`, :rcraw:`rightlabel.pad`, and :rcraw:`bottomlabel.pad`
        suptitle, figtitle : str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.
            This is an improvement on matplotlib's "super" title, which just
            centers the text between figure edges.
        suptitlepad : float, optional
            The padding between the super title and the axes content in arbitrary
            units (default is points). Default is :rcraw:`suptitle.pad`.

        Other parameters
        ----------------
        %(axes.other)s

        Important
        ---------
        The `abc`, `abcstyle`, `abcloc`, `titleloc`, and `titleabove` keywords and
        the various `pad` keywords are :ref:`configuration settings <ug_config>`.
        We explicitly document these arguments here because it is very common to change
        them. But many :ref:`other configuration settings <ug_format>` can be passed
        to ``format`` too.

        See also
        --------
        proplot.config.RcConfigurator.context
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        """
        # General figure settings
        # TODO: Work out awkward situation where "figure" settings applied on axes
        kw = rc.fill({'facecolor': 'figure.facecolor'}, context=True)
        self.figure.patch.update(kw)

        # General axes settings
        cycle = rc.get('axes.prop_cycle', context=True)
        if cycle is not None:
            self.set_prop_cycle(cycle)

        # Text positioning adjustments, applied at drawtime
        above = rc.get('title.above', context=True)
        if above is not None:
            self._title_above = above
        pad = rc.get('abc.titlepad', context=True)
        if pad is not None:
            self._abc_pad = pad
        pad = rc.get('title.pad', context=True)  # title
        if pad is not None:
            self._set_title_offset_trans(pad)
            self._title_pad = pad
        pad = rc.get('suptitle.pad', context=True)
        if pad is not None:
            self.figure._suptitle_pad = pad
        for side in tuple(self._label_pad):  # labels
            pad = rc.get(side + 'label.pad', context=True)
            if pad is not None:
                self._label_pad[side] = pad

        # Update a-b-c label and title(s)
        self._update_abc()
        self._update_title_all(
            title, l=ltitle, r=rtitle, c=ctitle,
            ul=ultitle, uc=uctitle, ur=urtitle,
            ll=lltitle, lc=lctitle, lr=lrtitle,
        )

        # Update 'super' labels
        suptitle = _not_none(figtitle=figtitle, suptitle=suptitle)
        rlabels = _not_none(rightlabels=rightlabels, rlabels=rlabels)
        blabels = _not_none(bottomlabels=bottomlabels, blabels=blabels)
        llabels = _not_none(rowlabels=rowlabels, leftlabels=leftlabels, llabels=llabels)
        tlabels = _not_none(collabels=collabels, toplabels=toplabels, tlabels=tlabels)
        self._update_super(
            suptitle, left=llabels, right=rlabels, top=tlabels, bottom=blabels
        )

    def area(self, *args, **kwargs):
        """
        Shorthand for `~matplotlib.axes.Axes.fill_between`.

        See also
        --------
        matplotlib.axes.Axes.fill_between
        proplot.axes.standardize_1d
        proplot.axes.fill_between_extras
        proplot.axes.apply_cycle
        """
        # NOTE: *Cannot* assign area = axes.Axes.fill_between because the
        # wrapper won't be applied and for some reason it messes up
        # automodsumm, which tries to put the matplotlib docstring on website
        return self.fill_between(*args, **kwargs)

    def areax(self, *args, **kwargs):
        """
        Shorthand for `~matplotlib.axes.Axes.fill_betweenx`.

        See also
        --------
        matplotlib.axes.Axes.fill_betweenx
        proplot.axes.standardize_1d
        proplot.axes.fill_betweenx_extras
        proplot.axes.apply_cycle
        """
        return self.fill_betweenx(*args, **kwargs)

    def boxes(self, *args, **kwargs):
        """
        Shorthand for `~matplotlib.axes.Axes.boxplot`.

        See also
        --------
        matplotlib.axes.Axes.boxplot
        proplot.axes.standardize_1d
        proplot.axes.boxplot_extras
        proplot.axes.apply_cycle
        """
        return self.boxplot(*args, **kwargs)

    def draw(self, renderer=None, *args, **kwargs):
        # Perform extra post-processing steps
        # NOTE: Used to have _reassign_title here (maybe _reassign_label too?)
        # but figured out it needs to get called by Figure spacing algorithm.
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
        bbox = super().get_tightbbox(renderer, *args, **kwargs)
        self._tight_bbox = bbox
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
        from .cartesian import CartesianAxes
        if not isinstance(self, CartesianAxes):
            warnings._warn_proplot(
                'Cannot adjust aspect ratio or ticks for non-Cartesian heatmap plot. '
                'Consider using pcolormesh() or pcolor() instead.'
            )
        else:
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
            transform = wrap._get_transform(self, transform)
        label = kwargs.pop('label', 'inset_axes')
        proj = _not_none(proj=proj, projection=projection)
        proj_kw = _not_none(proj_kw=proj_kw, projection_kw=projection_kw, default={})
        if basemap is not None:
            proj_kw['basemap'] = basemap

        # Inherit from current axes
        if proj is None:
            proj = self.name  # will have 'proplot_' prefix
            if proj_kw:
                warnings._warn_proplot(
                    'Inheriting projection from the main axes. '
                    f'Ignoring proj_kw keyword args: {proj_kw}'
                )
            if proj in ('proplot_cartopy', 'proplot_basemap'):
                m = copy.copy(self.projection)
                kwargs.setdefault('map_projection', m)

        # Create new projection
        elif proj == 'cartesian':
            proj = 'proplot_cartesian'
        elif proj == 'polar':
            proj = 'proplot_polar'
        else:
            m = constructor.Proj(proj, **proj_kw)
            kwargs.setdefault('map_projection', m)
            proj = 'proplot_' + m._proj_package

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
            The `zorder <https://matplotlib.org/stable/gallery/misc/zorder_demo.html>`__
            of the zoom lines. Should be greater than the zorder of
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

    def plotx(self, *args, **kwargs):
        """
        As with `~matplotlib.axes.Axes.plot` but interpret a single
        positional argument as *x* and multiple positional arguments
        as *y* and *x* (in that order).

        Parameters
        ----------
        *args, **kwargs
            Passed to `~matplotlib.axes.Axes.plot`.

        See also
        --------
        matplotlib.axes.Axes.plot
        proplot.axes.standardize_1d
        proplot.axes.indicate_error
        proplot.axes.apply_cycle
        """
        # NOTE: Arguments are standardized once we reach this block
        x, y, *args = args
        return super().plot(y, x, *args, **kwargs)

    @docstring.add_snippets
    def parametric(
        self, x, y, values=None, cmap=None, norm=None, *,
        interp=0, scalex=True, scaley=True, **kwargs
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
            The parametric coordinate can be indicated as a third positional
            argument or with the `values` or `levels` keywords.
        %(axes.cmap_norm)s
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
<https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line>`__.

        See also
        --------
        matplotlib.axes.Axes.plot
        proplot.axes.standardize_1d
        proplot.axes.apply_cmap
        """
        # Parse input
        # NOTE: Input *x* and *y* will have been standardized by _standardize_1d
        if values is None:
            raise ValueError('Values must be provided.')
        values = wrap._to_ndarray(values)
        ndim = tuple(_.ndim for _ in (x, y, values))
        size = tuple(_.size for _ in (x, y, values))
        if any(_ != 1 for _ in ndim):
            raise ValueError(f'Input coordinates must be 1D. Instead got dimensions {ndim}.')  # noqa: E501
        if any(_ != size[0] for _ in size):
            raise ValueError(f'Input coordinates must have identical size. Instead got sizes {size}.')  # noqa: E501

        # Interpolate values to allow for smooth gradations between values
        # (interp=False) or color switchover halfway between points
        # (interp=True). Then optionally interpolate the colormap values.
        # NOTE: The 'extras' wrapper handles input before ingestion by other wrapper
        # functions. *This* method is analogous to a native matplotlib method.
        if interp > 0:
            x_orig, y_orig, v_orig = x, y, values
            x, y, values = [], [], []
            for j in range(x_orig.shape[0] - 1):
                idx = slice(None)
                if j + 1 < x_orig.shape[0] - 1:
                    idx = slice(None, -1)
                x.extend(np.linspace(x_orig[j], x_orig[j + 1], interp + 2)[idx].flat)
                y.extend(np.linspace(y_orig[j], y_orig[j + 1], interp + 2)[idx].flat)
                values.extend(np.linspace(v_orig[j], v_orig[j + 1], interp + 2)[idx].flat)  # noqa: E501
            x, y, values = np.array(x), np.array(y), np.array(values)

        # Get coordinates and values for points to the 'left' and 'right' of joints
        coords = []
        levels = edges(values)
        for i in range(y.shape[0]):
            icoords = np.empty((3, 2))
            for j, arr in enumerate((x, y)):
                icoords[0, j] = arr[0] if i == 0 else 0.5 * (arr[i - 1] + arr[i])
                icoords[1, j] = arr[i]
                icoords[2, j] = arr[-1] if i + 1 == y.shape[0] else 0.5 * (arr[i + 1] + arr[i])  # noqa: E501
            coords.append(icoords)
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

    def scatterx(self, *args, **kwargs):
        """
        As with `~matplotlib.axes.Axes.scatter` but interpret a single
        positional argument as *x* and multiple positional arguments
        as *y* and *x* (in that order).

        Parameters
        ----------
        *args, **kwargs
            Passed to `~matplotlib.axes.Axes.scatter`.

        See also
        --------
        proplot.axes.standardize_1d
        matplotlib.axes.Axes.scatter
        proplot.axes.scatterx_extras
        """
        # NOTE: Arguments are standardized once we reach this block
        x, y, *args = args
        return super().scatter(y, x, *args, **kwargs)

    def violins(self, *args, **kwargs):
        """
        Shorthand for `~matplotlib.axes.Axes.violinplot`.

        See also
        --------
        matplotlib.axes.Axes.violinplot
        proplot.axes.standardize_1d
        proplot.axes.violinplot_extras
        proplot.axes.indicate_error
        proplot.axes.apply_cycle
        """
        return self.violinplot(*args, **kwargs)

    def _add_colorbar_legend(self, loc, obj, legend=False, **kwargs):
        """
        Queue up or replace objects for legends and list-of-artist style colorbars.
        """
        # Remove previous instances
        # NOTE: No good way to remove inset colorbars right now until the bounding
        # box and axes are merged into some colorbar subclass. Fine for now.
        d = self._legend_dict if legend else self._colorbar_dict
        if loc == 'fill':  # will be index in *parent* instead
            return
        if loc in d and not isinstance(d[loc], tuple):
            obj_prev = d.pop(loc)  # possibly pop a queued object
            if hasattr(self, 'legend_') and self.legend_ is obj_prev:
                self.legend_ = None  # was never added as artist
            elif legend:
                obj_prev.remove()  # remove legends and inner colorbars
        # Update queue or replace with instance
        if not isinstance(obj, tuple) or any(isinstance(_, mlegend.Legend) for _ in obj):  # noqa: E501
            d[loc] = obj
        else:
            handles, labels = obj
            handles_full, labels_full, kwargs_full = d.setdefault(loc, ([], [], {}))
            handles_full.extend(_not_none(handles, []))
            labels_full.extend(_not_none(labels, []))
            kwargs_full.update(kwargs)

    def _draw_colorbars_legends(self):
        """
        Draw the queued-up legends and colorbars. Wrapper funcs and legend func let
        user add handles to location lists with successive calls.
        """
        # WARNING: Passing empty list labels=[] to legend causes matplotlib
        # _parse_legend_args to search for everything. Ensure None if empty.
        for loc, colorbar in self._colorbar_dict.items():
            if not isinstance(colorbar, tuple):
                continue
            handles, labels, kwargs = colorbar
            self.colorbar(handles, labels or None, loc=loc, **kwargs)
        for loc, legend in self._legend_dict.items():
            if not isinstance(legend, tuple):
                continue
            elif any(isinstance(_, mlegend.Legend) for _ in legend):
                continue
            handles, labels, kwargs = legend
            self.legend(handles, labels or None, loc=loc, **kwargs)

    def _fill_colorbar_axes(self, length=None, **kwargs):
        """
        Return the axes and adjusted keyword args for a panel-filling colorbar.
        """
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

        # Draw colorbar axes within this one
        self._hide_panel()
        with self.figure._context_authorize_add_subplot():
            ax = self.figure.add_subplot(subplotspec, projection='proplot_cartesian')  # noqa: E501
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

        # Update default keyword args
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

        return ax, kwargs

    def _inset_colorbar_axes(
        self, loc=None, width=None, length=None, pad=None, frame=None, frameon=None,
        alpha=None, linewidth=None, edgecolor=None, facecolor=None, **kwargs
    ):
        """
        Return the axes and adjusted keyword args for an inset colorbar.
        """
        # Default properties
        cbwidth, cblength = width, length
        width, height = self.get_size_inches()
        frame = _not_none(frame=frame, frameon=frameon, default=rc['colorbar.frameon'])
        cbwidth = units(_not_none(cbwidth, rc['colorbar.insetwidth'])) / height
        cblength = units(_not_none(cblength, rc['colorbar.insetlength'])) / width
        extend = units(_not_none(kwargs.get('extendsize', None), rc['colorbar.insetextend']))  # noqa: E501
        pad = units(_not_none(pad, rc['colorbar.insetpad']))
        xpad, ypad = pad / width, pad / height

        # Get location in axes-relative coordinates
        # Bounds are x0, y0, width, height in axes-relative coordinates
        xspace = rc['xtick.major.size'] / 72
        if kwargs.get('label', None) or kwargs.get('title', None):
            xspace += 2.4 * rc['font.size'] / 72
        else:
            xspace += 1.2 * rc['font.size'] / 72
        xspace /= height  # space for labels
        if loc == 'upper right':
            ibounds = (1 - xpad - cblength, 1 - ypad - cbwidth)
            fbounds = (1 - 2 * xpad - cblength, 1 - 2 * ypad - cbwidth - xspace)
        elif loc == 'upper left':
            ibounds = (xpad, 1 - ypad - cbwidth)
            fbounds = (0, 1 - 2 * ypad - cbwidth - xspace)
        elif loc == 'lower left':
            ibounds = (xpad, ypad + xspace)
            fbounds = (0, 0)
        else:
            ibounds = (1 - xpad - cblength, ypad + xspace)
            fbounds = (1 - 2 * xpad - cblength, 0)
        ibounds = (*ibounds, cblength, cbwidth)  # inset axes
        fbounds = (*fbounds, 2 * xpad + cblength, 2 * ypad + cbwidth + xspace)  # frame

        # Make frame
        # NOTE: We do not allow shadow effects or fancy edges effect.
        # Also keep zorder same as with legend.
        if frame:
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
        locator = self._make_inset_locator(ibounds, self.transAxes)
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
            warnings._warn_proplot('Inset colorbars can only have ticks on the bottom.')
        kwargs.update({'orientation': 'horizontal', 'ticklocation': 'bottom'})
        kwargs.setdefault('maxn', 5)
        kwargs.setdefault('extendsize', extend)
        kwargs.update({'edgecolor': edgecolor, 'linewidth': linewidth})  # cbar edge

        return ax, kwargs

    def colorbar(
        self, mappable, values=None, *, loc=None, length=None, shrink=None, width=None,
        space=None, pad=None, queue=False, **kwargs
    ):
        """
        Add an *inset* colorbar or *outer* colorbar along the outside edge of
        the axes. See `~proplot.axes.colorbar_extras` for details.

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

        length : float or str, optional
            The colorbar length. For outer colorbars, units are relative to the
            axes width or height. Default is :rc:`colorbar.length`. For inset
            colorbars, units are interpreted by `~proplot.utils.units`. Default
            is :rc:`colorbar.insetlength`.
        shrink : float, optional
            Alias for `length`. This is included for consistency with
            `matplotlib.figure.Figure.colorbar`.
        width : float or str, optional
            The colorbar width. Units are interpreted by `~proplot.utils.units`.
            For outer colorbars, default is :rc:`colorbar.width`. For inset colorbars,
            default is :rc:`colorbar.insetwidth`.
        space : float or str, optional
            For outer colorbars only. The space between the colorbar and the
            main axes. Units are interpreted by `~proplot.utils.units`.
            When :rcraw:`tight` is ``True``, this is adjusted automatically.
            Otherwise, the default is :rc:`subplots.panelpad`.
        pad : float or str, optional
            For inset colorbars only. The space between the axes edge and the colorbar.
            Units are interpreted by `~proplot.utils.units`.
            Default is :rc:`colorbar.insetpad`.
        frame, frameon : bool, optional
            For inset colorbars only. Indicates whether to draw a "frame", just
            like `~matplotlib.axes.Axes.legend`. Default is :rc:`colorbar.frameon`.
        alpha, linewidth, edgecolor, facecolor : optional
            For inset colorbars only. Controls the transparency, edge width, edge color,
            and face color of the frame. Defaults are :rc:`colorbar.framealpha`,
            :rc:`axes.linewidth`, :rc:`axes.edgecolor`, and :rc:`axes.facecolor`.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~proplot.axes.colorbar_extras`.

        See also
        --------
        proplot.figure.Figure.colorbar
        proplot.axes.colorbar_extras
        """
        # TODO: Add option to pad the frame away from the axes edge
        # TODO: Get the 'best' inset colorbar location using the legend algorithm.
        length = _not_none(length=length, shrink=shrink)
        if loc != 'fill':
            loc = self._loc_translate(loc, 'colorbar')

        # Optionally add to queue
        if queue:
            obj = (mappable, values)
            kwargs.update({'width': width, 'length': length, 'space': space, 'pad': pad})  # noqa: E501
            return self._add_colorbar_legend(loc, obj, legend=False, **kwargs)

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            obj = ax.colorbar(mappable, values, loc='fill', length=length, **kwargs)
            self._add_colorbar_legend(loc, obj, legend=False)
            return obj

        # Generate colorbar axes
        if loc == 'fill':
            ax, kwargs = self._fill_colorbar_axes(length=length, **kwargs)
        else:
            ax, kwargs = self._inset_colorbar_axes(loc=loc, width=width, length=length, pad=pad, **kwargs)  # noqa: E501

        # Generate colorbar
        obj = wrap.colorbar_extras(ax, mappable, values, **kwargs)
        self._add_colorbar_legend(loc, obj, legend=False)  # possibly replace another
        return obj

    def legend(
        self, handles=None, labels=None, *,
        loc=None, width=None, space=None, queue=False, **kwargs
    ):
        """
        Add an *inset* legend or *outer* legend along the edge of the axes.
        See `~proplot.axes.legend_extras` for details.

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
        queue : bool, optional
            If ``True`` and `loc` is the same as an existing legend, the `handles`
            and `labels` are added to a queue and this function returns ``None``.
            This is used to "update" the same legend with successive ``ax.legend(...)``
            calls. If ``False`` (the default) and `loc` is the same as an existing
            legend, this function returns a `~matplotlib.legend.Legend` instance
            and the old legend is removed from the axes.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~proplot.axes.legend_extras`.

        See also
        --------
        proplot.figure.Figure.legend
        proplot.axes.legend_extras
        """
        if loc != 'fill':
            loc = self._loc_translate(loc, 'legend')
        if isinstance(loc, np.ndarray):
            loc = loc.tolist()

        # Optionally add to queue
        if queue:
            obj = (handles, labels)
            kwargs.update({'width': width, 'space': space})
            return self._add_colorbar_legend(loc, obj, legend=True, **kwargs)

        # Generate panel
        if loc in ('left', 'right', 'top', 'bottom'):
            ax = self.panel_axes(loc, width=width, space=space, filled=True)
            obj = ax.legend(handles, labels, loc='fill', **kwargs)
            self._add_colorbar_legend(loc, obj, legend=True)  # add to *this* axes
            return obj

        # Adjust settings
        if loc == 'fill':
            # Try to make handles and stuff flush against the axes edge
            self._hide_panel()
            kwargs.setdefault('borderaxespad', 0)
            frameon = _not_none(kwargs.get('frame'), kwargs.get('frameon'), rc['legend.frameon'])  # noqa: E501
            if not frameon:
                kwargs.setdefault('borderpad', 0)
            # Adjust location
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

        # Generate legend
        obj = wrap.legend_extras(self, handles, labels, loc=loc, **kwargs)
        self._add_colorbar_legend(loc, obj, legend=True)  # possibly replace another
        return obj

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
        for iax in (self, *(jax for side in panels for jax in self._panel_dict[side])):
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

    # Apply text wrapper
    text = wrap._apply_wrappers(
        maxes.Axes.text,
        wrap.text_extras,
    )

    # Apply 1D plotting command wrappers
    plot = wrap._apply_wrappers(
        maxes.Axes.plot,
        wrap.standardize_1d,
        wrap._plot_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    plotx = wrap._apply_wrappers(
        plotx,
        wrap.standardize_1d,
        wrap._plotx_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    step = wrap._apply_wrappers(
        maxes.Axes.step,
        wrap.standardize_1d,
        wrap.apply_cycle,
    )
    stem = wrap._apply_wrappers(
        maxes.Axes.stem,
        wrap.standardize_1d,
        wrap._stem_extras,
    )
    vlines = wrap._apply_wrappers(
        maxes.Axes.vlines,
        wrap.standardize_1d,
        wrap.vlines_extras,
        wrap.apply_cycle,
    )
    hlines = wrap._apply_wrappers(
        maxes.Axes.hlines,
        wrap.standardize_1d,
        wrap.hlines_extras,
        wrap.apply_cycle,
    )
    scatter = wrap._apply_wrappers(
        maxes.Axes.scatter,
        wrap.standardize_1d,
        wrap.scatter_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    scatterx = wrap._apply_wrappers(
        scatterx,
        wrap.standardize_1d,
        wrap.scatterx_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    bar = wrap._apply_wrappers(
        maxes.Axes.bar,
        wrap.standardize_1d,
        wrap.bar_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    barh = wrap._apply_wrappers(
        maxes.Axes.barh,
        wrap.standardize_1d,
        wrap.barh_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    hist = wrap._apply_wrappers(
        maxes.Axes.hist,
        wrap.standardize_1d,
        wrap.apply_cycle,
    )
    fill_between = wrap._apply_wrappers(
        maxes.Axes.fill_between,
        wrap.standardize_1d,
        wrap.fill_between_extras,
        wrap.apply_cycle,
    )
    fill_betweenx = wrap._apply_wrappers(
        maxes.Axes.fill_betweenx,
        wrap.standardize_1d,
        wrap.fill_betweenx_extras,
        wrap.apply_cycle,
    )
    boxplot = wrap._apply_wrappers(
        maxes.Axes.boxplot,
        wrap.standardize_1d,
        wrap.boxplot_extras,
        wrap.apply_cycle,
    )
    violinplot = wrap._apply_wrappers(
        maxes.Axes.violinplot,
        wrap.standardize_1d,
        wrap.violinplot_extras,
        wrap.indicate_error,
        wrap.apply_cycle,
    )
    pie = wrap._apply_wrappers(
        maxes.Axes.pie,
        wrap.standardize_1d,
        wrap.apply_cycle,
    )
    parametric = wrap._apply_wrappers(
        parametric,
        wrap.standardize_1d,
        wrap.apply_cmap,
    )
    hexbin = wrap._apply_wrappers(
        maxes.Axes.hexbin,
        wrap.standardize_1d,
        wrap.apply_cmap,
    )
    hist2d = wrap._apply_wrappers(
        maxes.Axes.hist2d,
        wrap.standardize_1d,
        wrap.apply_cmap,
    )

    # Apply 2D plotting command wrappers
    contour = wrap._apply_wrappers(
        maxes.Axes.contour,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    contourf = wrap._apply_wrappers(
        maxes.Axes.contourf,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    pcolor = wrap._apply_wrappers(
        maxes.Axes.pcolor,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    pcolormesh = wrap._apply_wrappers(
        maxes.Axes.pcolormesh,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    pcolorfast = wrap._apply_wrappers(
        maxes.Axes.pcolorfast,  # WARNING: not available in cartopy and basemap
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    streamplot = wrap._apply_wrappers(
        maxes.Axes.streamplot,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    quiver = wrap._apply_wrappers(
        maxes.Axes.quiver,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )
    barbs = wrap._apply_wrappers(
        maxes.Axes.barbs,
        wrap.standardize_2d,
        wrap.apply_cmap,
    )

    # Unstandardized commands
    tripcolor = wrap._apply_wrappers(
        maxes.Axes.tripcolor,
        wrap.apply_cmap,
    )
    tricontour = wrap._apply_wrappers(
        maxes.Axes.tricontour,
        wrap.apply_cmap,
    )
    tricontourf = wrap._apply_wrappers(
        maxes.Axes.tricontourf,
        wrap.apply_cmap,
    )
    spy = wrap._apply_wrappers(
        maxes.Axes.spy,
        wrap.apply_cmap,
    )
    imshow = wrap._apply_wrappers(
        maxes.Axes.imshow,
        wrap.apply_cmap,
    )
    matshow = wrap._apply_wrappers(
        maxes.Axes.matshow,
        wrap.apply_cmap,
    )
