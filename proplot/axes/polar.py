#!/usr/bin/env python3
"""
Polar axes using azimuth and radius instead of *x* and *y*.
"""
import inspect

import matplotlib.projections as mproj
import numpy as np

from .. import constructor
from .. import ticker as pticker
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import _not_none, _pop_rc, docstring, warnings
from . import plot, shared

__all__ = ['PolarAxes']


# Format docstring
_format_docstring = """
r0 : float, optional
    The radial origin. Default is ``0``.
theta0 : {'N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'}
    The zero azimuth location. Default is ``N``.
thetadir : {1, -1, 'anticlockwise', 'counterclockwise', 'clockwise'}, optional
    The positive azimuth direction. Clockwise corresponds to ``-1``
    and anticlockwise corresponds to ``1``. Default is ``1``.
thetamin, thetamax : float, optional
    The lower and upper azimuthal bounds in degrees. If
    ``thetamax != thetamin + 360``, this produces a sector plot.
thetalim : 2-tuple of float or None, optional
    Specifies `thetamin` and `thetamax` at once.
rmin, rmax : float, optional
    The inner and outer radial limits. If ``r0 != rmin``, this
    produces an annular plot.
rlim : 2-tuple of float or None, optional
    Specifies `rmin` and `rmax` at once.
rborder : bool, optional
    Whether to draw the polar axes border. Visibility of the "inner"
    radial spine and "start" and "end" azimuthal spines is controlled
    automatically by matplotlib.
thetagrid, rgrid, grid : bool, optional
    Whether to draw major gridlines for the azimuthal and radial axis.
    Use `grid` to toggle both.
thetagridminor, rgridminor, gridminor : bool, optional
    Whether to draw minor gridlines for the azimuthal and radial axis.
    Use `gridminor` to toggle both.
thetagridcolor, rgridcolor, gridcolor : color-spec, optional
    Color for the major and minor azimuthal and radial gridlines.
    Use `gridcolor` to set both at once.
thetalocator, rlocator : locator-spec, optional
    Used to determine the azimuthal and radial gridline positions.
    Passed to the `~proplot.constructor.Locator` constructor. Can be
    float, list of float, string, or `matplotlib.ticker.Locator` instance.
thetalines, rlines
    Aliases for `thetalocator`, `rlocator`.
thetalocator_kw, rlocator_kw : dict-like, optional
    The azimuthal and radial locator settings. Passed to
    `~proplot.constructor.Locator`.
thetaminorlocator, rminorlocator : optional
    As for `thetalocator`, `rlocator`, but for the minor gridlines.
thetaminorticks, rminorticks : optional
    Aliases for `thetaminorlocator`, `rminorlocator`.
thetaminorlocator_kw, rminorlocator_kw
    As for `thetalocator_kw`, `rlocator_kw`, but for the minor locator.
rlabelpos : float, optional
    The azimuth at which radial coordinates are labeled.
thetaformatter, rformatter : formatter-spec, optional
    Used to determine the azimuthal and radial label format.
    Passed to the `~proplot.constructor.Formatter` constructor.
    Can be string, list of string, or `matplotlib.ticker.Formatter`
    instance. Use ``[]``, ``'null'``, or ``'none'`` for no labels.
thetalabels, rlabels : optional
    Aliases for `thetaformatter`, `rformatter`.
thetaformatter_kw, rformatter_kw : dict-like, optional
    The azimuthal and radial label formatter settings. Passed to
    `~proplot.constructor.Formatter`.
labelpad : unit-spec, optional
    The padding between the axes edge and the radial and
    azimuthal labels. Default is :rcraw:`grid.labelpad`.
    %(units.pt)s
"""
docstring._snippet_manager['polar.format'] = _format_docstring


class PolarAxes(shared._SharedAxes, plot.PlotAxes, mproj.PolarAxes):
    """
    Axes subclass for plotting in polar coordinates. Adds the `~PolarAxes.format`
    method and overrides several existing methods.
    """
    _name = 'polar'

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        %(polar.format)s
        *args, **kwargs
            Passed to `proplot.axes.Axes`.

        Other parameters
        ----------------
        %(axes.format)s
        %(axes.rc)s

        See also
        --------
        PolarAxes.format
        proplot.axes.Axes
        proplot.axes.PlotAxes
        matplotlib.projections.PolarAxes
        proplot.figure.Figure.subplot
        proplot.figure.Figure.add_subplot
        """
        # Set tick length to zero so azimuthal labels are not too offset
        # Change default radial axis formatter but keep default theta one
        super().__init__(*args, **kwargs)
        self.yaxis.set_major_formatter(pticker.AutoFormatter())
        self.yaxis.isDefault_majfmt = True
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)

    def _update_formatter(self, x, *, formatter=None, formatter_kw=None):
        """
        Update the gridline label formatter.
        """
        # Tick formatter and toggling
        axis = getattr(self, x + 'axis')
        formatter_kw = formatter_kw or {}
        if formatter is not None:
            formatter = constructor.Formatter(formatter, **formatter_kw)  # noqa: E501
            axis.set_major_formatter(formatter)

    def _update_limits(self, x, *, min_=None, max_=None, lim=None):
        """
        Update the limits.
        """
        # Try to use public API where possible
        r = 'theta' if x == 'x' else 'r'
        if lim is not None:
            if min_ is not None or max_ is not None:
                warnings._warn_proplot(
                    f'Overriding {r}min={min_} and '
                    f'{r}max={max_} with {r}lim={lim}.'
                )
            min_, max_ = lim
        if min_ is not None:
            getattr(self, f'set_{r}min')(min_)
        if max_ is not None:
            getattr(self, f'set_{r}max')(max_)

    def _update_locators(
        self, x, *,
        locator=None, locator_kw=None, minorlocator=None, minorlocator_kw=None,
    ):
        """
        Update the gridline locator.
        """
        # TODO: Add minor tick 'toggling' as with cartesian axes?
        # NOTE: Must convert theta locator input to radians, then back to deg.
        r = 'theta' if x == 'x' else 'r'
        axis = getattr(self, x + 'axis')
        min_ = getattr(self, f'get_{r}min')()
        max_ = getattr(self, f'get_{r}max')()
        for i, (loc, loc_kw) in enumerate(
            zip((locator, minorlocator), (locator_kw, minorlocator_kw))
        ):
            if loc is None:
                continue
            # Get locator
            loc_kw = loc_kw or {}
            loc = constructor.Locator(loc, **loc_kw)
            # Sanitize values
            array = loc.tick_values(min_, max_)
            array = array[(array >= min_) & (array <= max_)]
            if x == 'x':
                array = np.deg2rad(array)
                if np.isclose(array[-1], min_ + 2 * np.pi):  # exclusive if 360 deg
                    array = array[:-1]
            # Assign fixed location
            loc = constructor.Locator(array)  # convert to FixedLocator
            if i == 0:
                axis.set_major_locator(loc)
            else:
                axis.set_minor_locator(loc)

    @docstring._snippet_manager
    def format(
        self, *, r0=None, theta0=None, thetadir=None,
        thetamin=None, thetamax=None, thetalim=None,
        rmin=None, rmax=None, rlim=None,
        thetagrid=None, rgrid=None,
        thetagridminor=None, rgridminor=None,
        thetagridcolor=None, rgridcolor=None,
        rlabelpos=None, rscale=None, rborder=None,
        thetalocator=None, rlocator=None, thetalines=None, rlines=None,
        thetalocator_kw=None, rlocator_kw=None,
        thetaminorlocator=None, rminorlocator=None, thetaminorlines=None, rminorlines=None,  # noqa: E501
        thetaminorlocator_kw=None, rminorlocator_kw=None,
        thetaformatter=None, rformatter=None, thetalabels=None, rlabels=None,
        thetaformatter_kw=None, rformatter_kw=None, labelpad=None,
        **kwargs
    ):
        """
        Modify axes limits, radial and azimuthal gridlines, and more. Note that
        all of the ``theta`` arguments are specified in degrees, not radians.

        Parameters
        ----------
        %(polar.format)s

        Other parameters
        ----------------
        %(axes.format)s
        %(figure.format)s
        %(axes.rc)s

        See also
        --------
        proplot.axes.Axes.format
        proplot.config.Configurator.context
        """
        # NOTE: Here we capture 'label.pad' rc argument normally used for
        # x and y axis labels as shorthand for 'tick.labelpad'.
        rc_kw, rc_mode = _pop_rc(kwargs)
        with rc.context(rc_kw, mode=rc_mode):
            # Not mutable default args
            thetalocator_kw = thetalocator_kw or {}
            thetaminorlocator_kw = thetaminorlocator_kw or {}
            thetaformatter_kw = thetaformatter_kw or {}
            rlocator_kw = rlocator_kw or {}
            rminorlocator_kw = rminorlocator_kw or {}
            rformatter_kw = rformatter_kw or {}

            # Flexible input
            thetalocator = _not_none(thetalines=thetalines, thetalocator=thetalocator)
            thetaformatter = _not_none(thetalabels=thetalabels, thetaformatter=thetaformatter)  # noqa: E501
            thetaminorlocator = _not_none(thetaminorlines=thetaminorlines, thetaminorlocator=thetaminorlocator)  # noqa: E501
            rlocator = _not_none(rlines=rlines, rlocator=rlocator)
            rformatter = _not_none(rlabels=rlabels, rformatter=rformatter)
            rminorlocator = _not_none(rminorlines=rminorlines, rminorlocator=rminorlocator)  # noqa: E501

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

            # Loop over axes
            for (
                x,
                grid,
                gridminor,
                gridcolor,
                min_, max_, lim,
                locator, locator_kw,
                formatter, formatter_kw,
                minorlocator, minorlocator_kw,
            ) in zip(
                ('x', 'y'),
                (thetagrid, rgrid),
                (thetagridminor, rgridminor),
                (thetagridcolor, rgridcolor),
                (thetamin, rmin), (thetamax, rmax), (thetalim, rlim),
                (thetalocator, rlocator), (thetalocator_kw, rlocator_kw),
                (thetaformatter, rformatter), (thetaformatter_kw, rformatter_kw),
                (thetaminorlocator, rminorlocator),
                (thetaminorlocator_kw, rminorlocator_kw),
            ):
                # Axis limits
                self._update_limits(x, min_=min_, max_=max_, lim=lim)

                # Axis tick settings
                # NOTE: Here use 'grid.labelpad' instead of 'tick.labelpad'. Default
                # offset for grid labels is larger than for tick labels.
                self._update_ticks(
                    x, grid=grid, gridminor=gridminor, gridcolor=gridcolor,
                    labelpad=labelpad, gridpad=True  # use grid.labelpad
                )

                # Axis locator
                self._update_locators(
                    x, locator=locator, locator_kw=locator_kw,
                    minorlocator=minorlocator, minorlocator_kw=minorlocator_kw
                )

                # Axis formatter
                self._update_formatter(
                    x, formatter=formatter, formatter_kw=formatter_kw
                )

        # Parent format method
        super().format(rc_kw=rc_kw, rc_mode=rc_mode, **kwargs)

    # Apply signature obfuscation after getting keys
    # NOTE: This is needed for __init__
    _format_signature = inspect.signature(format)
    format = docstring._obfuscate_kwargs(format)
