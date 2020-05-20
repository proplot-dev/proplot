#!/usr/bin/env python3
"""
The standard x-y axes used for most ProPlot figures.
"""
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from . import base
from .. import ticker as pticker
from .. import scale as pscale
from .. import constructor
from ..config import rc
from ..utils import units
from ..internals import ic  # noqa: F401
from ..internals import rcsetup, docstring, warnings, _not_none

__all__ = ['CartesianAxes']


_alt_doc = """
Return an axes in the same location as this one but whose {x} axis is on
the {x2}. This is an alias and more intuitive name for
`~CartesianAxes.twin{y}`, which generates two *{x}* axes with
a shared ("twin") *{y}* axes.

Also enforces the following settings:

* Places the old *{x}* axis on the {x1} and the new *{x}* axis
  on the {x2}.
* Makes the old {x2} spine invisible and the new {x1}, {y1},
  and {y2} spines invisible.
* Adjusts the *{x}* axis tick, tick label, and axis label positions
  according to the visible spine positions.
* Locks the old and new *{y}* axis limits and scales, and makes the new
  {y} axis labels invisible.

Parameters
----------
{xargs} : optional
    Passed to `Axes.format`.
{args} : optional
    Prepended with ``'{x}'`` and passed to `Axes.format`.
"""

_alt_kwargs = (  # TODO: More systematic approach?
    'label', 'locator', 'formatter', 'ticks', 'ticklabels',
    'minorlocator', 'minorticks', 'tickminor',
    'ticklen', 'tickrange', 'tickdir', 'ticklabeldir', 'tickrotation', 'wraprange',
    'bounds', 'margin', 'color', 'linewidth', 'grid', 'gridminor', 'gridcolor',
    'locator_kw', 'formatter_kw', 'minorlocator_kw', 'label_kw',
)

docstring.snippets['axes.altx'] = _alt_doc.format(
    x='x',
    x1='bottom',
    x2='top',
    y='y',
    y1='left',
    y2='right',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('x' + key for key in _alt_kwargs),
)

docstring.snippets['axes.alty'] = _alt_doc.format(
    x='y',
    x1='left',
    x2='right',
    y='x',
    y1='bottom',
    y2='top',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('y' + key for key in _alt_kwargs),
)

_dual_doc = """
Return a secondary *{x}* axis for denoting equivalent *{x}*
coordinates in *alternate units*.

Parameters
----------
arg : function, (function, function), or `~matplotlib.scale.ScaleBase`
    Used to transform units from the parent axis to the secondary axis.
    This can be a `~proplot.scale.FuncScale` itself or a function,
    (function, function) tuple, or `~matplotlib.scale.ScaleBase` instance used
    to *generate* a `~proplot.scale.FuncScale` (see
    `~proplot.scale.FuncScale` for details).
{args} : optional
    Prepended with ``'{x}'`` and passed to `Axes.format`.
"""

docstring.snippets['axes.dualx'] = _dual_doc.format(
    x='x',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('x' + key for key in _alt_kwargs),
)
docstring.snippets['axes.dualy'] = _dual_doc.format(
    x='y',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('y' + key for key in _alt_kwargs),
)

_twin_doc = """
Mimics the builtin `~matplotlib.axes.Axes.twin{y}` method.

Also enforces the following settings:

* Places the old *{x}* axis on the {x1} and the new *{x}* axis
  on the {x2}.
* Makes the old {x2} spine invisible and the new {x1}, {y1},
  and {y2} spines invisible.
* Adjusts the *{x}* axis tick, tick label, and axis label positions
  according to the visible spine positions.
* Locks the old and new *{y}* axis limits and scales, and makes the new
  {y} axis labels invisible.

Parameters
----------
{xargs} : optional
    Passed to `Axes.format`.
{args} : optional
    Prepended with ``'{x}'`` and passed to `Axes.format`.
"""

docstring.snippets['axes.twinx'] = _twin_doc.format(
    x='y', x1='left', x2='right',
    y='x', y1='bottom', y2='top',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('y' + key for key in _alt_kwargs),
)
docstring.snippets['axes.twiny'] = _twin_doc.format(
    x='x', x1='bottom', x2='top',
    y='y', y1='left', y2='right',
    args=', '.join(_alt_kwargs),
    xargs=', '.join('x' + key for key in _alt_kwargs),
)


def _parse_alt(x, kwargs):
    """
    Interpret keyword args passed to all "twin axis" methods so they
    can be passed to Axes.format.
    """
    kw_bad, kw_out = {}, {}
    for key, value in kwargs.items():
        if key in _alt_kwargs:
            kw_out[x + key] = value
        elif key[0] == x and key[1:] in _alt_kwargs:
            # NOTE: We permit both e.g. 'locator' and 'xlocator' because
            # while is more elegant and consistent with e.g. colorbar() syntax
            # but latter is more consistent and easier to use when refactoring.
            kw_out[key] = value
        elif key in rcsetup._rc_nodots:
            kw_out[key] = value
        else:
            kw_bad[key] = value
    if kw_bad:
        raise TypeError(f'Unexpected keyword argument(s): {kw_bad!r}')
    return kw_out


def _parse_rcloc(x, string):  # figures out string location
    """
    Convert the *boolean* "left", "right", "top", and "bottom" rc settings
    to a location string. Returns ``None`` if settings are unchanged.
    """
    if x == 'x':
        top = rc.get(f'{string}.top', context=True)
        bottom = rc.get(f'{string}.bottom', context=True)
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
    else:
        left = rc.get(f'{string}.left', context=True)
        right = rc.get(f'{string}.right', context=True)
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


class CartesianAxes(base.Axes):
    """
    Axes subclass for plotting in ordinary Cartesian coordinates.
    Adds the `~CartesianAxes.format` method and overrides several existing
    methods.
    """
    #: The registered projection name.
    name = 'cartesian'

    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        proplot.ui.subplots
        """
        # Impose default formatter
        super().__init__(*args, **kwargs)
        formatter = pticker.AutoFormatter()
        self.xaxis.set_major_formatter(formatter)
        self.yaxis.set_major_formatter(formatter)
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True
        # Custom attributes
        self._datex_rotated = False  # whether manual rotation has been applied
        self._dualy_arg = None  # for scaling units on opposite side of ax
        self._dualx_arg = None
        self._dualy_cache = None  # prevent excess _dualy_overrides calls
        self._dualx_cache = None

    def _altx_overrides(self):
        """
        Apply alternate *x* axis overrides.
        """
        # Unlike matplotlib API, we strong arm user into certain twin axes
        # settings... doesn't really make sense to have twin axes without this
        # NOTE: Could also use _panel_sharey_group = True to hide xaxis content
        # but instead we set entire axis to visible = False. Safer that way.
        if self._altx_child is not None:  # altx was called on this axes
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
        """
        Apply alternate *y* axis overrides.
        """
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

    def _apply_axis_sharing(self):
        """
        Enforce the "shared" axis labels and axis tick labels. If this is not called
        at drawtime, "shared" labels can be inadvertantly turned off.
        """
        # X axis
        axis = self.xaxis
        if self._sharex is not None:
            level = 3 if self._panel_sharex_group else self.figure._sharex
            if level > 0:
                axis.label.set_visible(False)
            if level > 2:
                # WARNING: Cannot set NullFormatter because shared axes share the
                # same axis.Ticker classes. Instead use the approach copied from
                # matplotlib subplots().
                axis.set_tick_params(which='both', labelbottom=False, labeltop=False)

        # Y axis
        axis = self.yaxis
        if self._sharey is not None:
            level = 3 if self._panel_sharey_group else self.figure._sharey
            if level > 0:
                axis.label.set_visible(False)
            if level > 2:
                axis.set_tick_params(which='both', labelleft=False, labelright=False)
        axis.set_minor_formatter(mticker.NullFormatter())

    def _datex_rotate(self):
        """
        Apply default rotation to datetime axis coordinates.
        """
        # NOTE: Rotation is done *before* horizontal/vertical alignment,
        # cannot change alignment with set_tick_params. Must apply to text
        # objects. fig.autofmt_date calls subplots_adjust, so cannot use it.
        if (
            not isinstance(self.xaxis.converter, mdates.DateConverter)
            or self._datex_rotated
        ):
            return
        rotation = rc['formatter.timerotation']
        kw = {'rotation': rotation}
        if rotation not in (0, 90, -90):
            kw['ha'] = ('right' if rotation > 0 else 'left')
        for label in self.xaxis.get_ticklabels():
            label.update(kw)
        self._datex_rotated = True  # do not need to apply more than once

    def _dualx_overrides(self):
        """
        Lock the child "dual" *x* axis limits to the parent.
        """
        # NOTE: We set the scale using private API to bypass application of
        # set_default_locators_and_formatters: only_if_default=True is critical
        # to prevent overriding user settings! We also bypass autoscale_view
        # because we set limits manually, and bypass child.stale = True
        # because that is done in call to set_xlim() below.
        arg = self._dualx_arg
        if arg is None:
            return
        scale = self.xaxis._scale
        olim = self.get_xlim()
        if (scale, *olim) == self._dualx_cache:
            return
        child = self._altx_child
        funcscale = pscale.FuncScale(arg, invert=True, parent_scale=scale)
        child.xaxis._scale = funcscale
        child._update_transScale()
        funcscale.set_default_locators_and_formatters(
            child.xaxis, only_if_default=True
        )
        nlim = list(map(funcscale.functions[1], np.array(olim)))
        if np.sign(np.diff(olim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]  # if function flips limits, so will set_xlim!
        child.set_xlim(nlim, emit=False)
        self._dualx_cache = (scale, *olim)

    def _dualy_overrides(self):
        """
        Lock the child "dual" *y* axis limits to the parent.
        """
        arg = self._dualy_arg
        if arg is None:
            return
        scale = self.yaxis._scale
        olim = self.get_ylim()
        if (scale, *olim) == self._dualy_cache:
            return
        child = self._alty_child
        funcscale = pscale.FuncScale(arg, invert=True, parent_scale=scale)
        child.yaxis._scale = funcscale
        child._update_transScale()
        funcscale.set_default_locators_and_formatters(
            child.yaxis, only_if_default=True
        )
        nlim = list(map(funcscale.functions[1], np.array(olim)))
        if np.sign(np.diff(olim)) != np.sign(np.diff(nlim)):
            nlim = nlim[::-1]
        child.set_ylim(nlim, emit=False)
        self._dualy_cache = (scale, *olim)

    def _sharex_setup(self, sharex):
        """
        Configure shared axes accounting for panels. The input is the
        'parent' axes, from which this one will draw its properties.
        """
        # Share panel across different subplots
        super()._sharex_setup(sharex)
        # Get sharing level
        level = (
            3 if self._panel_sharex_group and self._is_panel_parent_or_child(sharex)
            else self.figure._sharex
        )
        if level not in range(4):
            raise ValueError(
                'Invalid sharing level sharex={value!r}. '
                'Axis sharing level can be 0 (share nothing), '
                '1 (hide axis labels), '
                '2 (share limits and hide axis labels), or '
                '3 (share limits and hide axis and tick labels).'
            )
        if sharex in (None, self) or not isinstance(sharex, CartesianAxes):
            return
        # Builtin sharing features
        if level > 0:
            self._sharex = sharex
            text = self.xaxis.label.get_text()
            if text:
                sharex.xaxis.label.set_text(text)
        if level > 1:
            self._shared_x_axes.join(self, sharex)
            self.xaxis.major = sharex.xaxis.major
            self.xaxis.minor = sharex.xaxis.minor

    def _sharey_setup(self, sharey):
        """
        Configure shared axes accounting for panels. The input is the
        'parent' axes, from which this one will draw its properties.
        """
        # Share panel across different subplots
        super()._sharey_setup(sharey)
        # Get sharing level
        level = (
            3 if self._panel_sharey_group and self._is_panel_parent_or_child(sharey)
            else self.figure._sharey
        )
        if level not in range(4):
            raise ValueError(
                'Invalid sharing level sharey={value!r}. '
                'Axis sharing level can be 0 (share nothing), '
                '1 (hide axis labels), '
                '2 (share limits and hide axis labels), or '
                '3 (share limits and hide axis and tick labels).'
            )
        if sharey in (None, self) or not isinstance(sharey, CartesianAxes):
            return
        # Builtin features
        if level > 0:
            self._sharey = sharey
            text = self.yaxis.label.get_text()
            if text:
                sharey.yaxis.label.set_text(text)
        if level > 1:
            self._shared_y_axes.join(self, sharey)
            self.yaxis.major = sharey.yaxis.major
            self.yaxis.minor = sharey.yaxis.minor

    def _update_axis_labels(self, x='x', **kwargs):
        """
        Apply axis labels to the relevant shared axis. If spanning labels are toggled
        this keeps the labels synced for all subplots in the same row or column. Label
        positions will be adjusted at draw-time with figure._align_axislabels.
        """
        if x not in 'xy':
            return

        # Get axes in 3 step process
        # 1. Walk to parent if it is a main axes
        # 2. Get spanning main axes in this row or column (ignore short panel edges)
        # 3. Walk to parent if it exists (may be a panel long edge)
        # NOTE: Critical to apply labels to *shared* axes attributes rather
        # than testing extents or we end up sharing labels with twin axes.
        ax = self
        if getattr(self.figure, '_share' + x) > 0:
            share = getattr(ax, '_share' + x) or ax
            if not share._panel_parent:
                ax = share  # shared main axes are only ever one level deep

        # Get spanning axes
        axs = [ax]
        if getattr(ax.figure, '_span' + x):
            side = getattr(self, x + 'axis').get_label_position()
            if side in ('left', 'bottom'):
                axs = ax._get_side_axes(side, panels=False)

        # Update axes with label
        for ax in axs:
            ax = getattr(ax, '_share' + x) or ax  # defer to panel
            axis = getattr(ax, x + 'axis')
            axis.label.update(kwargs)

    @docstring.add_snippets
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
        xtickminor=None, ytickminor=None,
        xticklabeldir=None, yticklabeldir=None,
        xtickrange=None, ytickrange=None,
        xwraprange=None, ywraprange=None,
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
        xlinewidth=None, ylinewidth=None,
        xgridcolor=None, ygridcolor=None,
        xticklen=None, yticklen=None,
        xlabel_kw=None, ylabel_kw=None,
        xscale_kw=None, yscale_kw=None,
        xlocator_kw=None, ylocator_kw=None,
        xformatter_kw=None, yformatter_kw=None,
        xminorlocator_kw=None, yminorlocator_kw=None,
        patch_kw=None,
        **kwargs
    ):
        """
        Modify the *x* and *y* axis labels, tick locations, tick labels,
        axis scales, spine settings, and more. Unknown keyword arguments
        are passed to `Axes.format` and
        `~proplot.config.RcConfigurator.context`.

        Parameters
        ----------
        aspect : {'auto', 'equal'}, optional
            The aspect ratio mode. See `~matplotlib.axes.Axes.set_aspect`
            for details.
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
            `~proplot.scale.Scale` constructor. For example,
            ``xscale='log'`` applies logarithmic scaling, and
            ``xscale=('cutoff', 0.5, 2)`` applies a custom
            `~proplot.scale.CutoffScale`.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.scale.Scale`.
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
            Use `grid` to toggle both.
        xgridminor, ygridminor : bool, optional
            Whether to draw minor gridlines for the *x* and *y* axis.
            Use `gridminor` to toggle both.
        xticklabeldir, yticklabeldir : {'in', 'out'}
            Whether to place *x* and *y* axis tick label text inside
            or outside the axes.
        xlocator, ylocator : locator spec, optional
            Used to determine the *x* and *y* axis tick mark positions. Passed
            to the `~proplot.constructor.Locator` constructor.  Can be float,
            list of float, string, or `matplotlib.ticker.Locator` instance.
        xticks, yticks : optional
            Aliases for `xlocator`, `ylocator`.
        xlocator_kw, ylocator_kw : dict-like, optional
            The *x* and *y* axis locator settings. Passed to
            `~proplot.constructor.Locator`.
        xminorlocator, yminorlocator : optional
            As for `xlocator`, `ylocator`, but for the minor ticks.
        xminorticks, yminorticks : optional
            Aliases for `xminorlocator`, `yminorlocator`.
        xminorlocator_kw, yminorlocator_kw
            As for `xlocator_kw`, `ylocator_kw`, but for the minor locator.
        xformatter, yformatter : formatter spec, optional
            Used to determine the *x* and *y* axis tick label string format.
            Passed to the `~proplot.constructor.Formatter` constructor.
            Can be string, list of strings, or `matplotlib.ticker.Formatter`
            instance. Use ``[]`` or ``'null'`` for no ticks.
        xticklabels, yticklabels : optional
            Aliases for `xformatter`, `yformatter`.
        xformatter_kw, yformatter_kw : dict-like, optional
            The *x* and *y* axis formatter settings. Passed to
            `~proplot.constructor.Formatter`.
        xrotation, yrotation : float, optional
            The rotation for *x* and *y* axis tick labels. Default is ``0``
            for normal axes, :rc:`formatter.timerotation` for time
            *x* axes.
        xtickrange, ytickrange : (float, float), optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``(-1, 1)`` with
            axis range ``(-5, 5)`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1. See
            `~proplot.ticker.AutoFormatter` for details.
        xwraprange, ywraprange : (float, float), optional
            The *x* and *y* axis data ranges with which major tick mark
            values are *wrapped*. For example, the wrap range ``(0, 3)``
            causes the values 0 through 9 to be formatted as 0, 1, 2,
            0, 1, 2, 0, 1, 2, 0. See `~proplot.ticker.AutoFormatter` for details.
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
        xlinewidth, ylinewidth : color-spec, optional
            Line width for the *x* and *y* axis spines and major ticks.
            Default is :rc:`linewidth`. Use e.g. ``ax.format(linewidth=2)``
            to set for both axes.
        xgridcolor, ygridcolor : color-spec, optional
            Color for the *x* and *y* axis major and minor gridlines.
            Default is :rc:`grid.color`. Use e.g. ``ax.format(gridcolor='r')``
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
        %(axes.patch_kw)s

        Other parameters
        ----------------
        %(axes.other)s

        See also
        --------
        proplot.config.RcConfigurator.context
        proplot.axes.Axes.format

        Note
        ----
        If you plot something with a `datetime64 \
<https://docs.scipy.org/doc/numpy/reference/arrays.datetime.html>`__,
        `pandas.Timestamp`, `pandas.DatetimeIndex`, `datetime.date`,
        `datetime.time`, or `datetime.datetime` array as the *x* or *y* axis
        coordinate, the axis ticks and tick labels will be automatically
        formatted as dates.
        """
        rc_kw, rc_mode, kwargs = self._parse_format(**kwargs)
        with rc.context(rc_kw, mode=rc_mode):
            # Background patch
            kw_face = rc.fill(
                {
                    'facecolor': 'axes.facecolor',
                    'alpha': 'axes.alpha'
                },
                context=True,
            )
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
            xmargin = _not_none(xmargin, rc.get('axes.xmargin', context=True))
            ymargin = _not_none(ymargin, rc.get('axes.ymargin', context=True))
            xtickdir = _not_none(xtickdir, rc.get('xtick.direction', context=True))
            ytickdir = _not_none(ytickdir, rc.get('ytick.direction', context=True))
            xformatter = _not_none(xformatter=xformatter, xticklabels=xticklabels)
            yformatter = _not_none(yformatter=yformatter, yticklabels=yticklabels)
            xlocator = _not_none(xlocator=xlocator, xticks=xticks)
            ylocator = _not_none(ylocator=ylocator, yticks=yticks)
            xtickminor = _not_none(
                xtickminor, rc.get('xtick.minor.visible', context=True)
            )
            ytickminor = _not_none(
                ytickminor, rc.get('ytick.minor.visible', context=True)
            )
            xminorlocator = _not_none(
                xminorlocator=xminorlocator, xminorticks=xminorticks,
            )
            yminorlocator = _not_none(
                yminorlocator=yminorlocator, yminorticks=yminorticks,
            )

            # Grid defaults are more complicated
            grid = rc.get('axes.grid', context=True)
            which = rc.get('axes.grid.which', context=True)
            if which is not None or grid is not None:  # if *one* was changed
                axis = rc['axes.grid.axis']  # always need this property
                if grid is None:
                    grid = rc['axes.grid']
                elif which is None:
                    which = rc['axes.grid.which']
                xgrid = _not_none(
                    xgrid, grid and axis in ('x', 'both')
                    and which in ('major', 'both')
                )
                ygrid = _not_none(
                    ygrid, grid and axis in ('y', 'both')
                    and which in ('major', 'both')
                )
                xgridminor = _not_none(
                    xgridminor, grid and axis in ('x', 'both')
                    and which in ('minor', 'both')
                )
                ygridminor = _not_none(
                    ygridminor, grid and axis in ('y', 'both')
                    and which in ('minor', 'both')
                )

            # Sensible defaults for spine, tick, tick label, and label locs
            # NOTE: Allow tick labels to be present without ticks! User may
            # want this sometimes! Same goes for spines!
            xspineloc = _not_none(xloc=xloc, xspineloc=xspineloc,)
            yspineloc = _not_none(yloc=yloc, yspineloc=yspineloc,)
            xtickloc = _not_none(xtickloc, xspineloc, _parse_rcloc('x', 'xtick'))
            ytickloc = _not_none(ytickloc, yspineloc, _parse_rcloc('y', 'ytick'))
            xspineloc = _not_none(xspineloc, _parse_rcloc('x', 'axes.spines'))
            yspineloc = _not_none(yspineloc, _parse_rcloc('y', 'axes.spines'))
            if xtickloc != 'both':
                xticklabelloc = _not_none(xticklabelloc, xtickloc)
                xlabelloc = _not_none(xlabelloc, xticklabelloc)
                if xlabelloc not in (None, 'bottom', 'top'):  # e.g. "both"
                    xlabelloc = 'bottom'
            if ytickloc != 'both':
                yticklabelloc = _not_none(yticklabelloc, ytickloc)
                ylabelloc = _not_none(ylabelloc, yticklabelloc)
                if ylabelloc not in (None, 'left', 'right'):
                    ylabelloc = 'left'

            # Begin loop
            for (
                x, axis,
                label, color,
                linewidth, gridcolor,
                ticklen,
                margin, bounds,
                tickloc, spineloc,
                ticklabelloc, labelloc,
                grid, gridminor,
                tickminor, minorlocator,
                lim, reverse, scale,
                locator, tickrange,
                wraprange,
                formatter, tickdir,
                ticklabeldir, rotation,
                label_kw, scale_kw,
                locator_kw, minorlocator_kw,
                formatter_kw
            ) in zip(
                ('x', 'y'), (self.xaxis, self.yaxis),
                (xlabel, ylabel), (xcolor, ycolor),
                (xlinewidth, ylinewidth), (xgridcolor, ygridcolor),
                (xticklen, yticklen),
                (xmargin, ymargin), (xbounds, ybounds),
                (xtickloc, ytickloc), (xspineloc, yspineloc),
                (xticklabelloc, yticklabelloc), (xlabelloc, ylabelloc),
                (xgrid, ygrid), (xgridminor, ygridminor),
                (xtickminor, ytickminor), (xminorlocator, yminorlocator),
                (xlim, ylim), (xreverse, yreverse), (xscale, yscale),
                (xlocator, ylocator), (xtickrange, ytickrange),
                (xwraprange, ywraprange),
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
                    scale = constructor.Scale(scale, **scale_kw)
                    getattr(self, 'set_' + x + 'scale')(scale)

                # Is this a date axis?
                # NOTE: Make sure to get this *after* lims set!
                # See: https://matplotlib.org/api/units_api.html
                # And: https://matplotlib.org/api/dates_api.html
                # Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axis.py # noqa
                # The axis_date() method just applies DateConverter
                date = isinstance(axis.converter, mdates.DateConverter)

                # Fix spines
                kw = rc.fill(
                    {
                        'color': 'axes.edgecolor',
                        'linewidth': 'axes.linewidth',
                    },
                    context=True,
                )
                if color is not None:
                    kw['color'] = color
                if linewidth is not None:
                    kw['linewidth'] = linewidth
                sides = ('bottom', 'top') if x == 'x' else ('left', 'right')
                spines = [self.spines[side] for side in sides]
                for spine, side in zip(spines, sides):
                    # Line properties. Override if we're settings spine bounds
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
                                    f'Invalid {x} spine location {spineloc!r}. '
                                    'Options are: '
                                    + ', '.join(map(
                                        repr, (*sides, 'both', 'neither')
                                    )) + '.'
                                )
                    # Apply spine bounds
                    if bounds is not None and spine.get_visible():
                        spine.set_bounds(*bounds)
                    spine.update(kw)

                # Get available spines, needed for setting tick locations
                spines = [
                    side for side, spine in zip(sides, spines)
                    if spine.get_visible()
                ]

                # Tick and grid settings for major and minor ticks separately
                # Override is just a "new default", but user can override this
                for which, igrid in zip(('major', 'minor'), (grid, gridminor)):
                    # Tick properties
                    # NOTE: This loads xtick.major.size, xtick.major.width,
                    # xtick.major.pad, xtick.major.bottom, and xtick.major.top
                    # For all the x/y major/minor tick types
                    kwticks = rc.category(x + 'tick.' + which, context=True)
                    if kwticks is None:
                        kwticks = {}
                    else:
                        kwticks.pop('visible', None)  # invalid setting
                    if ticklen is not None:
                        kwticks['size'] = units(ticklen, 'pt')
                        if which == 'minor':
                            kwticks['size'] *= rc['ticklenratio']

                    # Grid style and toggling
                    name = 'grid' if which == 'major' else 'gridminor'
                    if igrid is not None:
                        axis.grid(igrid, which=which)
                    kwgrid = rc.fill(
                        {
                            'grid_color': name + '.color',
                            'grid_alpha': name + '.alpha',
                            'grid_linewidth': name + '.linewidth',
                            'grid_linestyle': name + '.linestyle',
                        },
                        context=True,
                    )
                    if gridcolor is not None:  # override for specific x/y axes
                        kw['grid_color'] = gridcolor
                    axis.set_tick_params(which=which, **kwgrid, **kwticks)

                # Tick and ticklabel properties that apply equally for major/minor lines
                # Weird issue causes set_tick_params to reset/forget grid is turned
                # on if you access tick.gridOn directly, instead of passing through
                # tick_params. Since gridOn is undocumented feature, don't use it.
                # So calling _format_axes() a second time will remove the lines.
                # First determine tick sides, avoiding situation where we draw ticks
                # on top of invisible spine.
                kw = {}
                loc2sides = {
                    None: None,
                    'both': sides,
                    'none': (),
                    'neither': (),
                }
                if bounds is not None and tickloc not in sides:
                    tickloc = sides[0]  # override to just one side
                ticklocs = loc2sides.get(tickloc, (tickloc,))
                if ticklocs is not None:
                    kw.update({side: side in ticklocs for side in sides})
                kw.update({side: False for side in sides if side not in spines})

                # Tick label sides
                # Will override to make sure only appear where ticks are
                ticklabellocs = loc2sides.get(ticklabelloc, (ticklabelloc,))
                if ticklabellocs is not None:
                    kw.update(
                        {'label' + side: (side in ticklabellocs) for side in sides}
                    )
                kw.update(  # override
                    {
                        'label' + side: False for side in sides
                        if side not in spines
                        or (ticklocs is not None and side not in ticklocs)
                    }
                )  # override

                # The axis label side
                if labelloc is None:
                    if ticklocs is not None:
                        options = [
                            side for side in sides
                            if side in ticklocs and side in spines
                        ]
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
                kw = rc.fill(
                    {
                        'labelcolor': 'tick.labelcolor',  # new props
                        'labelsize': 'tick.labelsize',
                        'color': x + 'tick.color',
                    },
                    context=True,
                )
                if color:
                    kw['color'] = color
                    kw['labelcolor'] = color

                # Tick label direction and rotation
                if tickdir == 'in':  # ticklabels should be much closer
                    kw['pad'] = 1.0
                if ticklabeldir == 'in':  # put tick labels inside the plot
                    tickdir = 'in'
                    kw['pad'] = -1.0 * sum(
                        rc[f'{x}tick.{key}']
                        for key in ('major.size', 'major.pad', 'labelsize')
                    )
                if tickdir is not None:
                    kw['direction'] = tickdir
                axis.set_tick_params(which='both', **kw)

                # Settings that can't be controlled by set_tick_params
                # Also set rotation and alignment here
                kw = rc.fill(
                    {
                        'fontfamily': 'font.family',
                        'weight': 'tick.labelweight'
                    },
                    context=True,
                )
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
                kw = rc.fill(
                    {
                        'color': 'axes.labelcolor',
                        'weight': 'axes.labelweight',
                        'fontsize': 'axes.labelsize',
                        'fontfamily': 'font.family',
                    },
                    context=True,
                )
                if label is not None:
                    kw['text'] = label
                if color:
                    kw['color'] = color
                kw.update(label_kw)
                if kw:  # NOTE: initially keep spanning labels off
                    self._update_axis_labels(x, **kw)

                # Major and minor locator
                # NOTE: Parts of API (dualxy) rely on minor tick toggling
                # preserving the isDefault_minloc setting. In future should
                # override the default matplotlib API minorticks_on!
                # NOTE: Unlike matplotlib API when "turning on" minor ticks
                # we *always* use the scale default, thanks to scale classes
                # refactoring with _ScaleBase. See Axes.minorticks_on.
                if locator is not None:
                    locator = constructor.Locator(locator, **locator_kw)
                    axis.set_major_locator(locator)
                    if isinstance(locator, mticker.IndexLocator):
                        tickminor = False  # 'index' minor ticks make no sense
                if minorlocator in (True, False):
                    warnings._warn_proplot(
                        f'You passed {x}minorticks={minorlocator}, but this '
                        'argument is used to specify tick *locations*. If '
                        'you just want to *toggle* minor ticks on and off, '
                        f'please use {x}tickminor=True or {x}tickminor=False.'
                    )
                    minorlocator = None
                if tickminor or minorlocator:
                    isdefault = minorlocator is None
                    if isdefault:
                        minorlocator = getattr(
                            axis._scale, '_default_minor_locator', None
                        )
                        if not minorlocator:
                            minorlocator = constructor.Locator('minor')
                    else:
                        minorlocator = constructor.Locator(
                            minorlocator, **minorlocator_kw
                        )
                    axis.set_minor_locator(minorlocator)
                    axis.isDefault_minloc = isdefault
                elif tickminor is not None and not tickminor:
                    # NOTE: Generally if you *enable* minor ticks on a dual
                    # axis, want to allow FuncScale updates to change the
                    # minor tick locators. If you *disable* minor ticks, do
                    # not want FuncScale applications to turn them on. So we
                    # allow below to set isDefault_minloc to False.
                    axis.set_minor_locator(constructor.Locator('null'))

                # Major formatter
                # NOTE: The only reliable way to disable ticks labels and then
                # restore them is by messing with the *formatter*, rather than
                # setting labelleft=False, labelright=False, etc.
                if (
                    formatter is not None
                    or tickrange is not None
                    or wraprange is not None
                ):
                    # Tick range
                    if tickrange is not None or wraprange is not None:
                        if formatter not in (None, 'auto'):
                            warnings._warn_proplot(
                                'The tickrange and autorange features require '
                                'proplot.AutoFormatter formatter. Overriding '
                                'input formatter.'
                            )
                        formatter = 'auto'
                        if tickrange is not None:
                            formatter_kw.setdefault('tickrange', tickrange)
                        if wraprange is not None:
                            formatter_kw.setdefault('wraprange', wraprange)

                    # Set the formatter
                    # Note some formatters require 'locator' as keyword arg
                    if formatter in ('date', 'concise'):
                        locator = axis.get_major_locator()
                        formatter_kw.setdefault('locator', locator)
                    formatter = constructor.Formatter(
                        formatter, date=date, **formatter_kw
                    )
                    axis.set_major_formatter(formatter)

                # Ensure no out-of-bounds ticks; set_smart_bounds() can fail
                # * Using set_bounds did not work, so instead just turn
                #   locators into fixed version.
                # * Most locators take no arguments in __call__, and some do
                #   not have tick_values method, so we just call them.
                if (
                    bounds is not None
                    or fixticks
                    or isinstance(formatter, mticker.FixedFormatter)
                    or axis.get_scale() == 'cutoff'
                ):
                    if bounds is None:
                        bounds = getattr(self, 'get_' + x + 'lim')()
                    locator = constructor.Locator([
                        x for x in axis.get_major_locator()()
                        if bounds[0] <= x <= bounds[1]
                    ])
                    axis.set_major_locator(locator)
                    locator = constructor.Locator([
                        x for x in axis.get_minor_locator()()
                        if bounds[0] <= x <= bounds[1]
                    ])
                    axis.set_minor_locator(locator)

            # Call parent
            if aspect is not None:
                self.set_aspect(aspect)
            super().format(**kwargs)

    @docstring.add_snippets
    def altx(self, **kwargs):
        """
        %(axes.altx)s
        """
        # Cannot wrap twiny() because we want to use CartesianAxes, not
        # matplotlib Axes. Instead use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py  # noqa
        if self._altx_child or self._altx_parent:
            raise RuntimeError('No more than *two* twin axes are allowed.')
        with self.figure._context_authorize_add_subplot():
            ax = self._make_twin_axes(sharey=self, projection='cartesian')
        ax.set_autoscaley_on(self.get_autoscaley_on())
        ax.grid(False)
        self._altx_child = ax
        ax._altx_parent = self
        self._altx_overrides()
        ax._altx_overrides()
        self.add_child_axes(ax)  # to facilitate tight layout
        self.figure._axstack.remove(ax)  # or gets drawn twice!
        ax.format(**_parse_alt('x', kwargs))
        return ax

    @docstring.add_snippets
    def alty(self, **kwargs):
        """
        %(axes.alty)s
        """
        # Docstring is programatically assigned below
        if self._alty_child or self._alty_parent:
            raise RuntimeError('No more than *two* twin axes are allowed.')
        with self.figure._context_authorize_add_subplot():
            ax = self._make_twin_axes(sharex=self, projection='cartesian')
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.grid(False)
        self._alty_child = ax
        ax._alty_parent = self
        self._alty_overrides()
        ax._alty_overrides()
        self.add_child_axes(ax)  # to facilitate tight layout
        self.figure._axstack.remove(ax)  # or gets drawn twice!
        ax.format(**_parse_alt('y', kwargs))
        return ax

    @docstring.add_snippets
    def dualx(self, arg, **kwargs):
        """
        %(axes.dualx)s
        """
        # NOTE: Matplotlib 3.1 has a 'secondary axis' feature. For the time
        # being, our version is more robust (see FuncScale) and simpler, since
        # we do not create an entirely separate _SecondaryAxis class.
        ax = self.altx(**kwargs)
        self._dualx_arg = arg
        self._dualx_overrides()
        return ax

    @docstring.add_snippets
    def dualy(self, arg, **kwargs):
        """
        %(axes.dualy)s
        """
        ax = self.alty(**kwargs)
        self._dualy_arg = arg
        self._dualy_overrides()
        return ax

    def draw(self, renderer=None, *args, **kwargs):
        # Perform extra post-processing steps
        # NOTE: This mimics matplotlib API, which calls identical
        # post-processing steps in both draw() and get_tightbbox()
        self._altx_overrides()
        self._alty_overrides()
        self._dualx_overrides()
        self._dualy_overrides()
        self._datex_rotate()
        self._apply_axis_sharing()
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        super().draw(renderer, *args, **kwargs)

    def get_tightbbox(self, renderer, *args, **kwargs):
        # Perform extra post-processing steps
        self._altx_overrides()
        self._alty_overrides()
        self._dualx_overrides()
        self._dualy_overrides()
        self._datex_rotate()
        self._apply_axis_sharing()
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        return super().get_tightbbox(renderer, *args, **kwargs)

    @docstring.add_snippets
    def twinx(self):
        """
        %(axes.twinx)s
        """
        return self.alty()

    @docstring.add_snippets
    def twiny(self):
        """
        %(axes.twiny)s
        """
        return self.altx()
