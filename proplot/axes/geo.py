#!/usr/bin/env python3
"""
Axes filled with cartographic projections.
"""
import numpy as np
import matplotlib.axes as maxes
import matplotlib.text as mtext
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from . import base
from .plot import (
    _basemap_norecurse, _basemap_redirect, _cmap_changer, _cycle_changer,
    _default_crs, _default_latlon, _default_transform,
    _fill_between_wrapper, _fill_betweenx_wrapper,
    _indicate_error, _plot_wrapper,
    _scatter_wrapper, _standardize_1d, _standardize_2d,
    _text_wrapper,
)
from .. import crs as pcrs
from .. import constructor
from .. import ticker as pticker
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import docstring, warnings, _version, _version_cartopy, _not_none
from ..utils import arange
try:
    from cartopy.mpl.geoaxes import GeoAxes as GeoAxesCartopy
    import cartopy.feature as cfeature
    import cartopy.mpl.ticker as cticker
    import cartopy.crs as ccrs
except ModuleNotFoundError:
    GeoAxesBase = object
    cfeature = cticker = ccrs = None
try:
    import mpl_toolkits.basemap as mbasemap
except ModuleNotFoundError:
    mbasemap = None

__all__ = ['GeoAxes', 'BasemapAxes', 'CartopyAxes']


def _circle_path(N=100):
    """
    Return a circle `~matplotlib.path.Path` used as the outline for polar
    stereographic, azimuthal equidistant, Lambert conformal, and gnomonic
    projections. This was developed from `this cartopy example \
<https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`__.
    """
    theta = np.linspace(0, 2 * np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)


class _GeoAxis(object):
    """
    Dummy axis used with longitude and latitude locators and for storing
    view limits on longitude latitude coordinates.
    """
    # NOTE: This is modeled after _DummyAxis in matplotlib/ticker.py, but we just
    # needs a couple methods used by simple, linear locators like `MaxNLocator`,
    # `MultipleLocator`, and `AutoMinorLocator`.
    def __init__(self, axes):
        self.axes = axes
        self.major = mticker.Ticker()
        self.minor = mticker.Ticker()
        self.isDefault_majfmt = True
        self.isDefault_majloc = True
        self.isDefault_minloc = True
        self._interval = None

    def get_scale(self):
        return 'linear'

    def get_tick_space(self):
        # Just use the long-standing default of nbins=9
        return 9

    def get_major_formatter(self):
        return self.major.formatter

    def get_major_locator(self):
        return self.major.locator

    def get_minor_locator(self):
        return self.minor.locator

    def get_majorticklocs(self):
        return self.major.locator()

    def get_minorticklocs(self):
        return self.minor.locator()

    def set_major_formatter(self, formatter, default=False):
        # NOTE: Cartopy <0.18 requires formatter axis is yaxis
        self.major.formatter = formatter
        formatter.set_axis(self.axis)
        self.isDefault_majfmt = default

    def set_major_locator(self, locator, default=False):
        self.major.locator = locator
        if self.major.formatter:
            self.major.formatter._set_locator(locator)
        locator.set_axis(self)
        self.isDefault_majloc = default

    def set_minor_locator(self, locator, default=False):
        self.minor.locator = locator
        locator.set_axis(self)
        self.isDefault_majfmt = default

    def set_view_interval(self, vmin, vmax):
        self._interval = (vmin, vmax)


class _LonAxis(_GeoAxis):
    """
    Axis with default longitude locator.
    """
    # NOTE: Basemap accepts tick formatters with drawmeridians(fmt=Formatter())
    # Try to use cartopy formatter if cartopy installed. Otherwise use
    # default builtin basemap formatting.
    def __init__(self, axes):
        self.axis = axes.xaxis  # the *actual* axis associated with this dummy one
        super().__init__(axes)
        if cticker is not None:
            self.set_major_formatter(cticker.LongitudeFormatter(), default=True)
        self.set_major_locator(pticker.LongitudeLocator(), default=True)
        self.set_minor_locator(mticker.AutoMinorLocator(), default=True)

    def get_view_interval(self):
        # NOTE: Proplot tries to set its *own* view intervals to avoid dateline
        # weirdness, but if cartopy.autoextent is False the interval will be
        # unset, so we are forced to use get_extent().
        interval = self._interval
        if interval is None:
            extent = self.axes.get_extent()
            interval = extent[:2]  # longitude extents
        return interval


class _LatAxis(_GeoAxis):
    """
    Axis with default latitude locator.
    """
    def __init__(self, axes):
        self._latmax = 90
        self.axis = axes.yaxis  # the *actual* axis associated with this dummy one
        super().__init__(axes)
        if cticker is not None:
            self.set_major_formatter(cticker.LatitudeFormatter(), default=True)
        self.set_major_locator(pticker.LatitudeLocator(), default=True)
        self.set_minor_locator(mticker.AutoMinorLocator(), default=True)

    def get_view_interval(self):
        interval = self._interval
        if interval is None:
            extent = self.axes.get_extent()
            interval = extent[2:]  # latitudes
        return interval

    def get_latmax(self):
        return self._latmax

    def set_latmax(self, latmax):
        self._latmax = latmax


class GeoAxes(base.Axes):
    """
    Axes subclass for plotting on cartographic projections.
    Adds the `~GeoAxes.format` method and overrides several existing methods.
    Subclassed by `CartopyAxes` and `BasemapAxes`, which respectively use
    cartopy `~cartopy.crs.Projection` and basemap `~mpl_toolkits.basemap.Basemap`
    objects to calculate projection coordinates.
    """
    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        proplot.ui.subplots
        proplot.axes.CartopyAxes
        proplot.axes.BasemapAxes
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

    @docstring.add_snippets
    def format(
        self, *,
        lonlim=None, latlim=None, boundinglat=None, grid=None,
        lonlines=None, lonlocator=None,
        latlines=None, latlocator=None, latmax=None,
        lonlines_kw=None, lonlocator_kw=None,
        latlines_kw=None, latlocator_kw=None,
        lonformatter=None, latformatter=None,
        lonformatter_kw=None, latformatter_kw=None,
        rotate_labels=None,
        labels=None, latlabels=None, lonlabels=None,
        patch_kw=None, **kwargs,
    ):
        """
        Modify the meridian and parallel labels, longitude and latitude map
        limits, geographic features, and more. Unknown keyword arguments are
        passed to `Axes.format` and `~proplot.config.rc_configurator.context`.

        Parameters
        ----------
        lonlim, latlim : (float, float), optional
            The approximate longitude and latitude boundaries of the map,
            applied with `~cartopy.mpl.geoaxes.GeoAxes.set_extent`. For
            cartopy axes only.
        boundinglat : float, optional
            The edge latitude for the circle bounding North Pole and
            South Pole-centered projections.
        grid : bool, optional
            Toggles meridian and parallel gridlines on and off. Default is
            :rc:`geogrid`.
        lonlines, latlines : float, list of float, or `~matplotlib.ticker.Locator`, \
optional
            If float, indicates the *spacing* of meridian and parallel gridlines. If
            list of float, indicates the exact meridian and parallel gridlines to draw.
            Otherwise, should be a `~matplotlib.ticker.Locator` instance.
        lonlocator, latlocator : optional
            Aliases for `lonlines`, `latlines`.
        lonlines_kw, latlines_kw : dict, optional
            Keyword argument dictionaries passed to the
            `~cartopy.mpl.ticker.LongitudeFormatter` and
            `~cartopy.mpl.ticker.LatitudeFormatter` formatters (respectively) for
            cartopy axes, or passed to `~mpl_toolkits.basemap.Basemap.drawmeridians`
            and `~mpl_toolkits.basemap.Basemap.drawparallels` methods (respectively)
            for basemap axes.
        lonlocator_kw, latlocator_kw : optional
            Aliases for `lonlines_kw`, `latlines_kw`.
        lonformatter, latformatter : `~matplotlib.ticker.Formatter`, optional
            `~matplotlib.ticker.Formatter` instances used to style longitude
            and latitude tick labels. For cartopy axes only.
        lonformatter_kw, latformatter_kw : optional
            Keyword arguments passed to `~cartopy.mpl.ticker.LongitudeFormatter`
            and `~cartopy.mpl.ticker.LatitudeFormatter`. Ignored if `lonformatter`
            or `latformatter` was provided. For cartopy axes only.
        rotate_labels : bool, optional
            Whether to rotate longitude and latitude gridline labels. For cartopy
            axes only. Default is :rc:`geogrid.rotatelabels`.
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
            2. A string indicating the side names, e.g. ``'lr'`` or ``'bt'``.
            3. A boolean 2-tuple indicating whether to draw labels on the
               ``(left, right)`` sides for longitudes, or ``(bottom, top)``
               sides for latitudes.
            4. A boolean 4-tuple indicating whether to draw labels on the
               ``(left, right, bottom, top)`` sides, as in the basemap
               `~mpl_toolkits.basemap.Basemap.drawmeridians` and
               `~mpl_toolkits.basemap.Basemap.drawparallels` methods.

        land, ocean, coast, rivers, lakes, borders, innerborders : bool, \
optional
            Toggles various geographic features. These are actually the
            :rcraw:`land`, :rcraw:`ocean`, :rcraw:`coast`, :rcraw:`rivers`,
            :rcraw:`lakes`, :rcraw:`borders`, and :rcraw:`innerborders`
            settings passed to `~proplot.config.rc_configurator.context`.
            The style can be modified using additional `rc` settings. For example,
            to change :rcraw:`land.color`, use ``ax.format(landcolor='green')``,
            and to change :rcraw:`land.zorder`, use ``ax.format(landzorder=4)``.
        %(axes.patch_kw)s

        Other parameters
        ----------------
        %(axes.other)s

        See also
        --------
        proplot.axes.Axes.format
        proplot.config.rc_configurator.context
        """
        # Format axes
        rc_kw, rc_mode, kwargs = self._parse_format(**kwargs)
        with rc.context(rc_kw, mode=rc_mode):
            # Parse alternative keyword args
            # TODO: Why isn't default latmax 80 respected sometimes?
            lonlines = _not_none(lonlines=lonlines, lonlocator=lonlocator)
            lonlines_user = lonlines is not None  # if True do not use LongitudeLocator
            lonlines = _not_none(lonlines, rc.get('geogrid.lonstep', context=True))
            latlines = _not_none(latlines=latlines, latlocator=latlocator)
            latlines_user = latlines is not None
            latlines = _not_none(latlines, rc.get('geogrid.latstep', context=True))
            lonlines_kw = _not_none(
                lonlines_kw=lonlines_kw, lonlocator_kw=lonlocator_kw, default={},
            )
            latlines_kw = _not_none(
                latlines_kw=latlines_kw, latlocator_kw=latlocator_kw, default={},
            )
            rotate_labels = _not_none(
                rotate_labels, rc.get('geogrid.rotatelabels', context=True)
            )
            lonformatter_kw = lonformatter_kw or {}
            latformatter_kw = latformatter_kw or {}
            latmax = _not_none(latmax, rc.get('geogrid.latmax', context=True))
            labels = _not_none(labels, rc.get('geogrid.labels', context=True))
            grid = _not_none(grid, rc.get('geogrid', context=True))
            if labels:
                lonlabels = _not_none(lonlabels, 1)
                latlabels = _not_none(latlabels, 1)

            # Get longitude lines
            if lonlines is not None:
                if np.iterable(lonlines):
                    lonlines = list(lonlines)
                else:
                    lonlines = self._get_lonlines(
                        step=lonlines, user=lonlines_user, **lonlines_kw
                    )

            # Get latitude lines. If latmax is changed we always need to reset latlines
            if latlines is not None or latmax is not None:
                if latlines is None:
                    latlines_prev = self._latlines_values
                    latlines = _not_none(latlines_prev, rc['geogrid.latstep'])
                if np.iterable(latlines):
                    latlines = list(latlines)
                else:
                    latlines = self._get_latlines(
                        step=latlines, latmax=latmax, user=latlines_user, **latlines_kw
                    )

            # Length-4 boolean arrays of whether and where to toggle labels
            # Format is [left, right, bottom, top]
            lonarray = self._parse_labels(lonlabels, lon=True)
            latarray = self._parse_labels(latlabels, lon=False)

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
            self._format_apply(
                patch_kw=patch_kw,
                boundinglat=boundinglat, lonlim=lonlim, latlim=latlim,
                lonlines=lonlines, latlines=latlines,
                lonlines_kw=lonlines_kw, latlines_kw=latlines_kw,
                lonformatter=lonformatter, latformatter=latformatter,
                lonformatter_kw=lonformatter_kw, latformatter_kw=latformatter_kw,
                rotate_labels=rotate_labels,
                latmax=latmax, lonarray=lonarray, latarray=latarray,
            )
            super().format(**kwargs)

    def _get_latlines(self, step, latmax=None):
        """
        Get latitude lines every `step` degrees.
        """
        # Latitudes gridlines, draw from -latmax to latmax unless result
        # would be asymmetrical across equator
        # NOTE: Basemap axes redraw *meridians* if they detect latmax was
        # explicitly changed, so important not to overwrite 'latmax'
        # with default value! Just need it for this calculation, then when
        # drawparallels is called will use self._latmax
        latmax = _not_none(latmax, self._latmax, rc['geogrid.latmax'])
        if latmax % step == -latmax % step:
            latlines = arange(-latmax, latmax, step)
        else:
            latlines = arange(0, latmax, step)
            if latlines[-1] != latmax:
                latlines = np.append(latlines, latmax)
            latlines = np.append(-latlines[::-1], latlines[1:])
        return list(latlines)

    def _get_lonlines(self, step, lon0=None):
        """
        Get longitude lines every `step` degrees.
        """
        # Longitude gridlines, draw relative to projection prime meridian
        # NOTE: We always generate gridlines array on first format call
        # because rc setting will be not None
        lon0 = _not_none(lon0, self._get_lon_0())
        lonlines = arange(lon0 - 180, lon0 + 180, step)
        lonlines = lonlines.astype(np.float64)
        if lonlines[-1] % 360 > 0:
            # Make sure the label appears on *right*, not on
            # top of the leftmost label.
            lonlines[-1] -= 1e-10
        else:
            # Formatter formats label as 1e-10... so there is
            # simply no way to put label on right. Just shift this
            # location off the map edge so parallels still extend
            # all the way to the edge, but label disappears.
            lonlines[-1] += 1e-10
        return list(lonlines)

    @staticmethod
    def _parse_labels(labels, lon):
        """
        Convert labels argument to length-4 boolean array.
        """
        if labels is None:
            return None
        if isinstance(labels, str):
            array = [0] * 4
            for idx, char in zip([0, 1, 2, 3], 'lrbt'):
                if char in labels:
                    array[idx] = 1
        else:
            array = np.atleast_1d(labels)
        if len(array) == 1:
            array = [*array, 0]  # default is to label bottom or left
        if len(array) == 2:
            if lon:
                array = [0, 0, *array]
            else:
                array = [*array, 0, 0]
        elif len(array) != 4:
            name = 'lon' if lon else 'lat'
            raise ValueError(f'Invalid {name}label spec: {labels}.')
        return array


class CartopyAxes(GeoAxes, GeoAxesCartopy):
    """
    Axes subclass for plotting
    `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__ projections.
    Makes ``transform=cartopy.crs.PlateCarree()`` the default for all plotting methods,
    enforces `global extent <https://stackoverflow.com/a/48956844/4970632>`__
    for most projections by default,  and draws `circular boundaries \
<https://scitools.org.uk/cartopy/docs/latest/gallery/always_circular_stereo.html>`__
    around polar azimuthal, stereographic, and Gnomonic projections bounded at
    the equator by default.
    """
    #: The registered projection name.
    name = 'cartopy'
    _circle_points = 100  # number of points for drawing circle map boundary

    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~cartopy.crs.Projection`
            The `~cartopy.crs.Projection` instance.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~cartopy.mpl.geoaxes.GeoAxes`.

        See also
        --------
        proplot.axes.GeoAxes
        proplot.constructor.Proj
        """
        # GeoAxes initialization. Note that critical attributes like
        # outline_patch needed by _format_apply are added before it is called.
        import cartopy.crs as ccrs
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError(
                'GeoAxes requires map_projection=cartopy.crs.Projection.'
            )
        super().__init__(*args, map_projection=map_projection, **kwargs)

        # Zero out ticks so gridlines are not offset
        for axis in (self.xaxis, self.yaxis):
            axis.set_tick_params(which='both', size=0)

        # Set extent and boundary extent for projections
        # The default bounding latitude is set in _format_apply
        # NOTE: set_global does not mess up non-global projections like OSNI
        if hasattr(self, 'set_boundary') and isinstance(
            self.projection, (
                ccrs.NorthPolarStereo, ccrs.SouthPolarStereo,
                pcrs.NorthPolarGnomonic, pcrs.SouthPolarGnomonic,
                pcrs.NorthPolarAzimuthalEquidistant,
                pcrs.NorthPolarLambertAzimuthalEqualArea,
                pcrs.SouthPolarAzimuthalEquidistant,
                pcrs.SouthPolarLambertAzimuthalEqualArea
            )
        ):
            path = _circle_path(self._circle_points)
            self.set_boundary(path, transform=self.transAxes)
        else:
            self.set_global()

    def _get_lon_0(self):
        """Get the central longitude."""
        return self.projection.proj4_params.get('lon_0', 0)

    def _get_lonlines(self, step, user=False, **kwargs):
        """Get longitude locator given the input step."""
        if not user and _version_cartopy >= _version('0.18'):
            from cartopy.mpl import ticker
            return ticker.LongitudeLocator(**kwargs)
        else:
            return super()._get_lonlines(step)

    def _get_latlines(self, step, latmax=None, user=False, **kwargs):
        """Get latitude locator given the input step."""
        # NOTE: After cartopy v0.18 meridian lines are no longer limited by the
        # most southerly and northerly parallel lines, for both the automatic
        # latitude locator and FixedLocator.
        if not user and _version_cartopy >= _version('0.18'):
            from cartopy.mpl import ticker
            return ticker.LatitudeLocator(**kwargs)
        else:
            return super()._get_latlines(step, latmax=latmax)

    def _format_apply(
        self, *, patch_kw,
        lonlim, latlim, boundinglat,
        lonlines, latlines, lonlines_kw, latlines_kw,
        lonformatter, latformatter, lonformatter_kw, latformatter_kw,
        rotate_labels,
        latmax, lonarray, latarray,
    ):
        """
        Apply formatting to cartopy axes. Extra kwargs are used to update proj4 params.
        """
        latmax  # prevent U100 error (cartopy handles 'latmax' automatically)
        lonlines_kw, latlines_kw  # preven U100 error (these were already applied)
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs
        from cartopy.mpl import ticker

        # Gridliner labels names
        def _toggle_labels(gl, left, right, bottom, top):
            if _version_cartopy >= _version('0.18'):  # cartopy >= 0.18
                left_labels = 'left_labels'
                right_labels = 'right_labels'
                bottom_labels = 'bottom_labels'
                top_labels = 'top_labels'
            else:  # cartopy < 0.18
                left_labels = 'ylabels_left'
                right_labels = 'ylabels_right'
                bottom_labels = 'xlabels_bottom'
                top_labels = 'xlabels_top'
            if left is not None:
                setattr(gl, left_labels, left)
            if right is not None:
                setattr(gl, right_labels, right)
            if bottom is not None:
                setattr(gl, bottom_labels, bottom)
            if top is not None:
                setattr(gl, top_labels, top)

        # Cartopy 0.18 monkey patch. This fixes issue where we get overlapping
        # gridlines on dateline. See the "nx -= 1" line in Gridliner._draw_gridliner
        # TODO: Submit cartopy PR. This is awful but necessary for quite a while if
        # the time between v0.17 and v0.18 is any indication.
        def _draw_gridliner(self, *args, **kwargs):
            result = type(self)._draw_gridliner(self, *args, **kwargs)
            if _version_cartopy == _version('0.18'):
                lon_lim, _ = self._axes_domain()
                if abs(np.diff(lon_lim)) == abs(np.diff(self.crs.x_limits)):
                    for collection in self.xline_artists:
                        if not getattr(collection, '_cartopy_fix', False):
                            collection.get_paths().pop(-1)
                            collection._cartopy_fix = True
            return result

        # Cartopy < 0.18 monkey patch. This is part of filtering valid label coordinates
        # to values between lon_0 - 180 and lon_0 + 180.
        def _axes_domain(self, *args, **kwargs):
            x_range, y_range = type(self)._axes_domain(self, *args, **kwargs)
            if _version_cartopy < _version('0.18'):
                lon_0 = self.axes.projection.proj4_params.get('lon_0', 0)
                x_range = np.asarray(x_range) + lon_0
            return x_range, y_range

        # Cartopy < 0.18 gridliner method monkey patch. Always print number in range
        # (180W, 180E). We choose #4 of the following choices 3 choices (see Issue #78):
        # 1. lonlines go from -180 to 180, but get double 180 labels at dateline
        # 2. lonlines go from -180 to e.g. 150, but no lines from 150 to dateline
        # 3. lonlines go from lon_0 - 180 to lon_0 + 180 mod 360, but results
        #    in non-monotonic array causing double gridlines east of dateline
        # 4. lonlines go from lon_0 - 180 to lon_0 + 180 monotonic, but prevents
        #    labels from being drawn outside of range (-180, 180)
        def _add_gridline_label(self, value, axis, upper_end):
            if _version_cartopy < _version('0.18'):
                if axis == 'x':
                    value = (value + 180) % 360 - 180
            return type(self)._add_gridline_label(self, value, axis, upper_end)

        # Initial gridliner object which ProPlot passively modifies
        # NOTE: The 'xpadding' and 'ypadding' props were introduced in v0.16
        # with default 5 points, then set to default None in v0.18.
        # TODO: Cartopy has had two formatters for a while but we use newer one
        # https://github.com/SciTools/cartopy/pull/1066
        if not self._gridliners:
            gl = self.gridlines(crs=ccrs.PlateCarree())
            gl._draw_gridliner = _draw_gridliner.__get__(gl)  # apply monkey patch
            gl._axes_domain = _axes_domain.__get__(gl)
            gl._add_gridline_label = _add_gridline_label.__get__(gl)
            gl.xlines = gl.ylines = False
            gl.xformatter = ticker.LongitudeFormatter()
            gl.yformatter = ticker.LatitudeFormatter()
            if _version_cartopy < _version('0.18'):  # necessary... for some reason...
                gl.xformatter.axis = self.xaxis
                gl.yformatter.axis = self.yaxis
            _toggle_labels(gl, False, False, False, False)
        gl = self._gridliners[0]

        # Projection extent
        # NOTE: They may add this as part of set_xlim and set_ylim in future
        # See: https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638  # noqa
        # WARNING: The set_extent method tries to set a *rectangle* between
        # the *4* (x,y) coordinate pairs (each corner), so something like
        # (-180, 180, -90, 90) will result in *line*, causing error!
        proj = self.projection.proj4_params['proj']
        north = isinstance(self.projection, (
            ccrs.NorthPolarStereo, pcrs.NorthPolarGnomonic,
            pcrs.NorthPolarAzimuthalEquidistant,
            pcrs.NorthPolarLambertAzimuthalEqualArea
        ))
        south = isinstance(self.projection, (
            ccrs.SouthPolarStereo, pcrs.SouthPolarGnomonic,
            pcrs.SouthPolarAzimuthalEquidistant,
            pcrs.SouthPolarLambertAzimuthalEqualArea
        ))
        if north or south:
            if lonlim is not None or latlim is not None:
                warnings._warn_proplot(
                    f'{proj!r} extent is controlled by "boundinglat", '
                    f'ignoring lonlim={lonlim!r} and latlim={latlim!r}.'
                )
            if self._boundinglat is None:
                if isinstance(self.projection, pcrs.NorthPolarGnomonic):
                    boundinglat = 30
                elif isinstance(self.projection, pcrs.SouthPolarGnomonic):
                    boundinglat = -30
                else:
                    boundinglat = 0
            if boundinglat is not None and boundinglat != self._boundinglat:
                eps = 1e-10  # bug with full -180, 180 range when lon_0 != 0
                lat0 = (90 if north else -90)
                lon0 = self.projection.proj4_params.get('lon_0', 0)
                extent = [lon0 - 180 + eps, lon0 + 180 - eps, boundinglat, lat0]
                self.set_extent(extent, crs=ccrs.PlateCarree())
                self._boundinglat = boundinglat
        else:
            if boundinglat is not None:
                warnings._warn_proplot(
                    f'{proj!r} extent is controlled by "lonlim" and "latlim", '
                    f'ignoring boundinglat={boundinglat!r}.'
                )
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

        # Gridline collection properties including axes.axisbelow-mimicking property
        kw = rc.fill({
            'alpha': 'geogrid.alpha',
            'color': 'geogrid.color',
            'linewidth': 'geogrid.linewidth',
            'linestyle': 'geogrid.linestyle',
        }, context=True)
        axisbelow = rc.get('geogrid.axisbelow', context=True)
        if axisbelow is not None:
            if axisbelow is True:
                zorder = 0.5
            elif axisbelow is False:
                zorder = 2.5
            elif axisbelow == 'lines':
                zorder = 1.5
            else:
                raise ValueError(f'Unexpected geogrid.axisbelow value {axisbelow!r}.')
            kw['zorder'] = zorder
        gl.collection_kwargs.update(kw)

        # Special gridline properties
        pad = rc.get('geogrid.labelpad', context=True)
        if pad is not None:
            gl.xpadding = gl.ypadding = pad
        loninline = rc.get('geogrid.loninline', context=True)
        if loninline is not None:
            gl.x_inline = loninline
        latinline = rc.get('geogrid.latinline', context=True)
        if latinline is not None:
            gl.y_inline = latinline

        # Gridline longitudes and latitudes
        eps = 1e-10
        if lonlines is not None:
            if isinstance(lonlines, mticker.Locator):
                gl.xlines = True
                gl.xlocator = lonlines
            else:
                # As of 0.18 lines no longer have to be in [lon_0 - 180, lon_0 + 180]
                if _version_cartopy >= _version('0.18'):
                    lonlines = (np.asarray(lonlines) + 180) % 360 - 180
                gl.xlines = bool(len(lonlines))
                if gl.xlines:
                    gl.xlocator = mticker.FixedLocator(lonlines)
        if latlines is not None:
            if isinstance(latlines, mticker.Locator):
                gl.ylines = True
                gl.ylocator = latlines
            else:
                gl.ylines = bool(len(latlines))
                if gl.ylines:
                    if latlines[0] == -90:
                        latlines[0] += eps
                    if latlines[-1] == 90:
                        latlines[-1] -= eps
                    gl.ylocator = mticker.FixedLocator(latlines)

        # Gridline label format
        if rotate_labels is not None:
            gl.rotate_labels = rotate_labels  # ignored in cartopy <0.18
        if lonformatter:
            gl.xformatter = lonformatter
        elif lonformatter_kw:
            gl.xformatter = ticker.LongitudeFormatter(**lonformatter_kw)
        if latformatter:
            gl.yformatter = latformatter
        elif latformatter_kw:
            gl.yformatter = ticker.LatitudeFormatter(**latformatter_kw)

        # Gridline label toggling
        # Issue warning instead of error!
        if _version_cartopy < _version('0.18'):
            if not isinstance(self.projection, (ccrs.Mercator, ccrs.PlateCarree)):
                if latarray is not None and any(latarray):
                    warnings._warn_proplot(
                        'Cannot add gridline labels to cartopy '
                        f'{type(self.projection).__name__} projection.'
                    )
                    latarray = [0] * 4
                if lonarray is not None and any(lonarray):
                    warnings._warn_proplot(
                        'Cannot add gridline labels to cartopy '
                        f'{type(self.projection).__name__} projection.'
                    )
                    lonarray = [0] * 4
        latarray = latarray or (None,) * 4
        lonarray = lonarray or (None,) * 4
        _toggle_labels(gl, *latarray[:2], *lonarray[2:])

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
        reso = rc['reso']
        reso = constructor.CARTOPY_RESOS.get(reso, None)
        if reso is None:
            raise ValueError(
                f'Invalid resolution {reso!r}. Options are: '
                + ', '.join(map(repr, constructor.CARTOPY_RESOS)) + '.'
            )
        for name, args in constructor.CARTPOY_FEATURES.items():
            # Get feature
            if not rc[name]:  # toggled
                continue
            if getattr(self, '_' + name, None):  # already drawn
                continue
            feat = cfeature.NaturalEarthFeature(*args, reso)
            # For 'lines', need to specify edgecolor and facecolor
            # See: https://github.com/SciTools/cartopy/issues/803
            kw = rc.category(name)  # do not omit uncached props
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
            'facecolor': 'geoaxes.facecolor',
            'alpha': 'geoaxes.facealpha',
        }, context=True)
        kw_edge = rc.fill({
            'edgecolor': 'geoaxes.edgecolor',
            'linewidth': 'geoaxes.linewidth',
        }, context=True)
        kw_face.update(patch_kw or {})
        self.background_patch.update(kw_face)
        self.outline_patch.update(kw_edge)

    def _apply_axis_sharing(self):
        """
        No-op for now. In future this will hide meridian and parallel
        labels for rectangular projections with axis sharing.
        """
        pass

    def get_tightbbox(self, renderer, *args, **kwargs):
        # Perform extra post-processing steps
        # For now this just draws the gridliners
        self._apply_axis_sharing()
        if self.get_autoscale_on() and self.ignore_existing_data_limits:
            self.autoscale_view()
        if (
            getattr(self.background_patch, 'reclip', None)
            and hasattr(self.background_patch, 'orig_path')
        ):
            clipped_path = self.background_patch.orig_path.clip_to_bbox(self.viewLim)
            self.background_patch._path = clipped_path

        # Adjust location
        if _version_cartopy >= _version('0.18'):
            self.patch._adjust_location()

        # Apply aspect
        self.apply_aspect()
        for gl in self._gridliners:
            if _version_cartopy >= _version('0.18'):
                gl._draw_gridliner(renderer=renderer)
            else:
                gl._draw_gridliner(background_patch=self.background_patch)

        # Remove gridliners
        if _version_cartopy < _version('0.18'):
            self._gridliners = []

        return super().get_tightbbox(renderer, *args, **kwargs)

    @property
    def projection(self):
        """
        The `~cartopy.crs.Projection` instance associated with this axes.
        """
        return self._map_projection

    @projection.setter
    def projection(self, map_projection):
        import cartopy.crs as ccrs
        if not isinstance(map_projection, ccrs.CRS):
            raise ValueError('Projection must be a cartopy.crs.CRS instance.')
        self._map_projection = map_projection

    # Wrapped methods
    # TODO: Remove this duplication!
    if GeoAxesCartopy is not object:
        text = _text_wrapper(
            GeoAxesCartopy.text
        )
        plot = _default_transform(_plot_wrapper(_standardize_1d(
            _indicate_error(_cycle_changer(GeoAxesCartopy.plot))
        )))
        scatter = _default_transform(_scatter_wrapper(_standardize_1d(
            _indicate_error(_cycle_changer(GeoAxesCartopy.scatter))
        )))
        fill_between = _fill_between_wrapper(_standardize_1d(_cycle_changer(
            GeoAxesCartopy.fill_between
        )))
        fill_betweenx = _fill_betweenx_wrapper(_standardize_1d(_cycle_changer(
            GeoAxesCartopy.fill_betweenx
        )))
        contour = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.contour
        )))
        contourf = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.contourf
        )))
        pcolor = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.pcolor
        )))
        pcolormesh = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.pcolormesh
        )))
        quiver = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.quiver
        )))
        streamplot = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.streamplot
        )))
        barbs = _default_transform(_standardize_2d(_cmap_changer(
            GeoAxesCartopy.barbs
        )))
        tripcolor = _default_transform(_cmap_changer(
            GeoAxesCartopy.tripcolor
        ))
        tricontour = _default_transform(_cmap_changer(
            GeoAxesCartopy.tricontour
        ))
        tricontourf = _default_transform(_cmap_changer(
            GeoAxesCartopy.tricontourf
        ))
        get_extent = _default_crs(
            GeoAxesCartopy.get_extent
        )
        set_extent = _default_crs(
            GeoAxesCartopy.set_extent
        )
        set_xticks = _default_crs(
            GeoAxesCartopy.set_xticks
        )
        set_yticks = _default_crs(
            GeoAxesCartopy.set_yticks
        )


class BasemapAxes(GeoAxes):
    """
    Axes subclass for plotting `~mpl_toolkits.basemap` projections. The
    `~mpl_toolkits.basemap.Basemap` instance is added as the
    `~BasemapAxes.projection` attribute, but you do not have to work with it
    directly -- plotting methods like `matplotlib.axes.Axes.plot` and
    `matplotlib.axes.Axes.contour` are redirected to the corresponding methods on
    the `~mpl_toolkits.basemap.Basemap` instance. Also ``latlon=True`` is passed
    to plotting methods by default.
    """
    #: The registered projection name.
    name = 'basemap'
    _proj_non_rectangular = (  # do not use axes spines as boundaries
        'ortho', 'geos', 'nsper',
        'moll', 'hammer', 'robin',
        'eck4', 'kav7', 'mbtfpq',
        'sinu', 'vandg',
        'npstere', 'spstere', 'nplaea',
        'splaea', 'npaeqd', 'spaeqd',
    )

    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `Axes`.

        See also
        --------
        proplot.axes.GeoAxes
        proplot.constructor.Proj
        """
        # WARNING: Investigated whether Basemap.__init__() could be called
        # twice with updated proj kwargs to modify map bounds after creation
        # and python immmediately crashes. Do not try again.
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
                'BasemapAxes requires map_projection=basemap.Basemap'
            )
        self._map_projection = map_projection
        self._map_boundary = None
        self._has_recurred = False  # use this to override plotting methods
        super().__init__(*args, **kwargs)

    def _get_lon_0(self, step=5):
        """Get the central longitude."""
        return step * round(self.projection.lonmin / step) + 180

    def _get_lonlines(self, step, user=False, **kwargs):
        """Get longitude line locations given the input step."""
        user, kwargs  # prevent U100 error (this is used in cartopy subclass)
        # Locations do not have to wrap around like they do in cartopy
        lonlines = super()._get_lonlines(step)
        lonlines = lonlines[:-1]
        return lonlines

    def _get_latlines(self, step, latmax=None, user=False, **kwargs):
        """Get latitude line locations given the input step."""
        user, kwargs  # prevent U100 error (these are used in cartopy subclass)
        return super()._get_latlines(step, latmax=latmax)

    def _format_apply(
        self, *, patch_kw,
        lonlim, latlim, boundinglat,
        lonlines, latlines, lonlines_kw, latlines_kw,
        lonformatter, latformatter, lonformatter_kw, latformatter_kw,
        rotate_labels,
        latmax, lonarray, latarray,
    ):
        """
        Apply changes to the basemap axes. Extra kwargs are used
        to update the proj4 params.
        """
        # Ignored arguments
        rotate_labels  # avoid U100 error (this arg is simply ignored)
        if lonformatter or lonformatter_kw:
            warnings._warn_proplot(
                f'Ignoring formatter arguments lonformatter={lonformatter} '
                f'with lonformatter_kw={lonformatter_kw} for basemap axes.'
            )
        if latformatter or latformatter_kw:
            warnings._warn_proplot(
                f'Ignoring formatter arguments latformatter={latformatter} '
                f'with latformatter_kw={latformatter_kw} for basemap axes.'
            )
        if (
            lonlim is not None
            or latlim is not None
            or boundinglat is not None
        ):
            warnings._warn_proplot(
                f'Got lonlim={lonlim!r}, latlim={latlim!r}, '
                f'boundinglat={boundinglat!r}, but you cannot "zoom into" a '
                'basemap projection after creating it. Pass proj_kw in your '
                'call to subplots with any of the following basemap keywords: '
                "'boundinglat', 'llcrnrlon', 'llcrnrlat', "
                "'urcrnrlon', 'urcrnrlat', 'llcrnrx', 'llcrnry', "
                "'urcrnrx', 'urcrnry', 'width', or 'height'."
            )

        # Map boundary
        # * First have to *manually replace* the old boundary by just
        #   deleting the original one
        # * If boundary is drawn successfully should be able to call
        #   self.projection._mapboundarydrawn.set_visible(False) and
        #   edges/fill color disappear
        # * For now will enforce that map plots *always* have background
        #   whereas axes plots can have transparent background
        kw_face = rc.fill({
            'facecolor': 'geoaxes.facecolor',
            'alpha': 'geoaxes.facealpha',
        }, context=True)
        kw_edge = rc.fill({
            'linewidth': 'geoaxes.linewidth',
            'edgecolor': 'geoaxes.edgecolor',
        }, context=True)
        kw_face.update(patch_kw or {})
        self.axesPatch = self.patch  # bugfix or something
        if self.projection.projection in self._proj_non_rectangular:
            self.patch.set_alpha(0)  # make patch invisible
            if not self.projection._mapboundarydrawn:
                # set fill_color to 'none' to make transparent
                p = self.projection.drawmapboundary(ax=self)
            else:
                p = self.projection._mapboundarydrawn
            p.update(kw_face)
            p.update(kw_edge)
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
        })
        tkw = rc.fill({
            'color': 'geogrid.color',
            'fontsize': 'geogrid.labelsize',
        })
        if lonarray is not None:  # change from lrbt to lrtb
            lonarray[2:] = lonarray[2:][::-1]
        if latarray is not None:  # change from lrbt to lrtb
            latarray[2:] = latarray[2:][::-1]

        # Parallel lines
        ilatmax = _not_none(latmax, self._latmax)
        if latlines is not None or latmax is not None or latarray is not None:
            for obj in self._iter_lines(self._latlines):
                obj.set_visible(False)
            latlines = _not_none(latlines, self._latlines_values)
            latarray = _not_none(latarray, self._latlines_labels, [0] * 4)
            p = self.projection.drawparallels(
                latlines, latmax=ilatmax, labels=latarray, ax=self, **latlines_kw,
            )
            self._latlines = p
            for obj in self._iter_lines(p):
                # Tried passing clip_on to the below, but it does nothing
                # Must set for lines created after the fact
                if isinstance(obj, mtext.Text):
                    obj.update(tkw)
                else:
                    obj.update(lkw)

        # Meridian lines
        # NOTE: Although it is not stated, latmax only affects drawmeridians for
        # north polar and south polar projections! Ignored otherwise!
        if lonlines is not None or latmax is not None or lonarray is not None:
            for obj in self._iter_lines(self._lonlines):
                obj.set_visible(False)
            lonlines = _not_none(lonlines, self._lonlines_values)
            lonarray = _not_none(lonarray, self._lonlines_labels, [0] * 4)
            p = self.projection.drawmeridians(
                lonlines, latmax=ilatmax, labels=lonarray, ax=self, **lonlines_kw,
            )
            self._lonlines = p
            for obj in self._iter_lines(p):
                if isinstance(obj, mtext.Text):
                    obj.update(tkw)
                else:
                    obj.update(lkw)

        # Geography
        # TODO: Allow setting the zorder.
        # NOTE: Also notable are drawcounties, blumarble, drawlsmask,
        # shadedrelief, and etopo methods.
        for name, method in constructor.BASEMAP_FEATURES.items():
            if not rc[name]:  # toggled
                continue
            if getattr(self, f'_{name}', None):  # already drawn
                continue
            kw = rc.category(name)
            feat = getattr(self.projection, method)(ax=self)
            if isinstance(feat, (list, tuple)):  # list of artists?
                for obj in feat:
                    obj.update(kw)
            else:
                feat.update(kw)
            setattr(self, '_' + name, feat)

    @staticmethod
    def _iter_lines(dict_):
        """Iterate over longitude latitude lines."""
        dict_ = dict_ or {}
        for pi in dict_.values():
            for pj in pi:
                for obj in pj:
                    yield obj

    @property
    def projection(self):
        """
        The `~mpl_toolkits.basemap.Basemap` instance associated with this axes.
        """
        return self._map_projection

    @projection.setter
    def projection(self, map_projection):
        import mpl_toolkits.basemap as mbasemap
        if not isinstance(map_projection, mbasemap.Basemap):
            raise ValueError('Projection must be a basemap.Basemap instance.')
        self._map_projection = map_projection

    # Wrapped methods
    plot = _basemap_norecurse(_default_latlon(_plot_wrapper(_standardize_1d(
        _indicate_error(_cycle_changer(_basemap_redirect(maxes.Axes.plot)))
    ))))
    scatter = _basemap_norecurse(_default_latlon(_scatter_wrapper(_standardize_1d(
        _indicate_error(_cycle_changer(_basemap_redirect(maxes.Axes.scatter)))
    ))))
    contour = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.contour)
    ))))
    contourf = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.contourf)
    ))))
    pcolor = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.pcolor)
    ))))
    pcolormesh = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.pcolormesh)
    ))))
    quiver = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.quiver)
    ))))
    streamplot = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.streamplot)
    ))))
    barbs = _basemap_norecurse(_default_latlon(_standardize_2d(_cmap_changer(
        _basemap_redirect(maxes.Axes.barbs)
    ))))
    hexbin = _basemap_norecurse(_standardize_1d(_cmap_changer(
        _basemap_redirect(maxes.Axes.hexbin)
    )))
    imshow = _basemap_norecurse(_cmap_changer(
        _basemap_redirect(maxes.Axes.imshow)
    ))
