#!/usr/bin/env python
"""
An axes used to jointly format Cartesian and polar axes.
"""
# NOTE: We could define these in base.py but idea is projection-specific formatters
# should never be defined on the base class. Might add to this class later anyway.
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import _pop_kwargs
from ..utils import _fontsize_to_pt, _not_none, units


class _SharedAxes(object):
    """
    Mix-in class with methods shared between `~proplot.axes.CartesianAxes`
    and `~proplot.axes.PolarAxes`.
    """
    def _update_background(self, x=None, **kwargs):
        """
        Update the background patch and spines.
        """
        # Update the background patch
        kw_face, kw_edge = self._get_background_props(**kwargs)
        self.patch.update(kw_face)
        if x is None:
            opts = self.spines
        elif x == 'x':
            opts = ('bottom', 'top', 'inner', 'polar')
        else:
            opts = ('left', 'right', 'start', 'end')
        for opt in opts:
            self.spines.get(opt, {}).update(kw_edge)

        # Update the tick colors
        axis = 'both' if x is None else x
        edgecolor = kw_edge.pop('edgecolor', None)
        if edgecolor is not None:
            self.tick_params(axis=axis, which='both', color=edgecolor)

        # Update the tick widths
        # TODO: Either exclude case where 'linewidth' was retrieved from
        # 'axes.linewidth' or make tick width a child of that setting?
        # TODO: The same logic is used inside config.py to scale tick widths
        # by tick ratios and to zero-out tick length. Share with helper func?
        linewidth = kw_edge.pop('linewidth', None)
        if linewidth is not None:
            kw = {'length': 0} if linewidth == 0 else {}
            self.tick_params(axis=axis, which='major', width=linewidth, **kw)
            ratio = rc['tick.widthratio']
            self.tick_params(axis=axis, which='minor', width=linewidth * ratio, **kw)

    def _update_ticks(
        self, x, *, grid=None, gridminor=None, gridcolor=None, gridpad=None,
        ticklen=None, tickdir=None, ticklabeldir=None, tickcolor=None, labelpad=None
    ):
        """
        Update the gridlines and labels. Set `gridpad` to ``True`` to use grid padding.
        """
        # Apply tick settings with tick_params when possible
        x = _not_none(x, 'x')
        kwtext = self._get_ticklabel_props(x)
        kwextra = _pop_kwargs(kwtext, 'weight', 'family')
        kwtext = {'label' + key: value for key, value in kwtext.items()}
        for b, which in zip((grid, gridminor), ('major', 'minor')):
            # Tick properties
            # NOTE: Must make 'tickcolor' overwrite 'labelcolor' or else 'color'
            # passed to __init__ will not apply correctly. Annoying but unavoidable
            kwticks = self._get_tick_props(x, which=which)
            if labelpad is not None:
                kwticks['pad'] = labelpad
            if tickcolor is not None:
                kwticks['color'] = tickcolor
                if rc._context_mode == 2:
                    kwtext.setdefault('labelcolor', tickcolor)
                else:
                    kwtext['labelcolor'] = tickcolor
            if ticklen is not None:
                kwticks['size'] = units(ticklen, 'pt')
                if which == 'minor':
                    kwticks['size'] *= rc['tick.lenratio']
            if gridpad:  # use grid.labelpad instead of tick.labelpad
                kwticks.pop('pad', None)
                pad = rc.find('grid.labelpad', context=True)
                if pad is not None:
                    kwticks['pad'] = units(pad, 'pt')

            # Gridline properties
            # NOTE: Internally ax.grid() passes gridOn to ax.tick_params() but this
            # is undocumented and might have weird side effects. Just use ax.grid()
            b = self._get_gridline_toggle(b, axis=x, which=which)
            if b is not None:
                self.grid(b, axis=x, which=which)
            kwlines = self._get_gridline_props(native=True, which=which)
            if gridcolor is not None:
                kwlines['grid_color'] = gridcolor

            # Apply tick and gridline properties
            self.tick_params(axis=x, which=which, **kwticks, **kwlines, **kwtext)

        # Tick and tick label direction with padding corrections
        # NOTE: The 'tick label direction' is right now just a cartesian thing
        kwdir = {}
        if tickdir == 'in':  # ticklabels should be much closer
            kwdir['pad'] = 1.0
        if ticklabeldir == 'in':  # put tick labels inside the plot
            tickdir = 'in'
            kwdir['pad'] = (
                - rc[f'{x}tick.major.size']
                - rc[f'{x}tick.major.pad']
                - _fontsize_to_pt(rc[f'{x}tick.labelsize'])
            )
        if tickdir is not None:
            kwdir['direction'] = tickdir
        self.tick_params(axis=x, which='both', **kwdir)

        # Apply settings that can't be controlled with tick_params
        if kwextra:
            axis = getattr(self, x + 'axis')
            for obj in axis.get_ticklabels():
                obj.update(kwextra)
