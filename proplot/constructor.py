#!/usr/bin/env python3
"""
The constructor functions used to build class instances from simple shorthand arguments.
"""
# NOTE: These functions used to be in separate files like crs.py and
# ticker.py but makes more sense to group them together to ensure usage is
# consistent and so online documentation is easier to understand. Also in
# future version classes will not be imported into top-level namespace. This
# change will be easier to do with all constructor functions in separate file.
# NOTE: Used to include the raw variable names that define string keys as
# part of documentation, but this is redundant and pollutes the namespace.
# User should just inspect docstrings, use trial-error, or see online tables.
import os
import re
from functools import partial
from numbers import Number

import cycler
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.projections.polar as mpolar
import matplotlib.scale as mscale
import matplotlib.ticker as mticker
import numpy as np

from . import colors as pcolors
from . import crs as pcrs
from . import scale as pscale
from . import ticker as pticker
from .config import rc
from .internals import ic  # noqa: F401
from .internals import (
    _not_none,
    _pop_props,
    _version,
    _version_cartopy,
    _version_mpl,
    warnings,
)
from .utils import to_rgba

try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    Basemap = object
try:
    import cartopy.crs as ccrs
    import cartopy.mpl.ticker as cticker
    from cartopy.crs import CRS
except ModuleNotFoundError:
    CRS = ccrs = cticker = object

__all__ = [
    'Colormap', 'Colors', 'Cycle', 'Norm',
    'Formatter', 'Locator', 'Scale', 'Proj',
]

# Dictionary of possible normalizers. See `Norm` for a table.
NORMS = {
    'none': mcolors.NoNorm,
    'null': mcolors.NoNorm,
    'div': pcolors.DivergingNorm,
    'diverging': pcolors.DivergingNorm,
    'segmented': pcolors.LinearSegmentedNorm,
    'log': mcolors.LogNorm,
    'linear': mcolors.Normalize,
    'power': mcolors.PowerNorm,
    'symlog': mcolors.SymLogNorm,
    'zero': pcolors.DivergingNorm,  # deprecated
    'midpoint': pcolors.DivergingNorm,  # deprecated
    'segments': pcolors.LinearSegmentedNorm,  # deprecated
}
if hasattr(mcolors, 'TwoSlopeNorm'):
    NORMS['twoslope'] = mcolors.TwoSlopeNorm

# Mapping of strings to `~matplotlib.ticker.Locator` classes. See
# `Locator` for a table."""
LOCATORS = {
    'none': mticker.NullLocator,
    'null': mticker.NullLocator,
    'auto': mticker.AutoLocator,
    'log': mticker.LogLocator,
    'maxn': mticker.MaxNLocator,
    'linear': mticker.LinearLocator,
    'multiple': mticker.MultipleLocator,
    'fixed': mticker.FixedLocator,
    'index': mticker.IndexLocator,
    'symlog': mticker.SymmetricalLogLocator,
    'logit': mticker.LogitLocator,
    'minor': mticker.AutoMinorLocator,
    'date': mdates.AutoDateLocator,
    'microsecond': mdates.MicrosecondLocator,
    'second': mdates.SecondLocator,
    'minute': mdates.MinuteLocator,
    'hour': mdates.HourLocator,
    'day': mdates.DayLocator,
    'weekday': mdates.WeekdayLocator,
    'month': mdates.MonthLocator,
    'year': mdates.YearLocator,
    'lon': partial(pticker.LongitudeLocator, dms=False),
    'lat': partial(pticker.LatitudeLocator, dms=False),
    'deglon': partial(pticker.LongitudeLocator, dms=False),
    'deglat': partial(pticker.LatitudeLocator, dms=False),
}
if _version_cartopy >= _version('0.18'):
    # NOTE: This only makes sense when paired with degree-minute-second formatter
    # NOTE: We copied cartopy locators because they are short and necessary
    # for determining both cartopy and basemap tick locations. We did *not* copy
    # formatter because they are long and we have nice, simpler alternatives of
    # deglon and deglat.
    LOCATORS['dms'] = partial(pticker._DegreeLocator, dms=True)
    LOCATORS['dmslon'] = partial(pticker.LongitudeLocator, dms=True)
    LOCATORS['dmslat'] = partial(pticker.LatitudeLocator, dms=True)
if hasattr(mpolar, 'ThetaLocator'):
    LOCATORS['theta'] = mpolar.ThetaLocator

# Mapping of strings to `~matplotlib.ticker.Formatter` classes. See
# `Formatter` for a table.
# NOTE: Critical to use SimpleFormatter for cardinal formatters rather than
# AutoFormatter because latter fails with Basemap formatting.
# NOTE: Define cartopy longitude/latitude formatters with dms=True because that
# is their distinguishing feature relative to proplot formatter.
FORMATTERS = {  # note default LogFormatter uses ugly e+00 notation
    'auto': pticker.AutoFormatter,
    'frac': pticker.FracFormatter,
    'sci': pticker.SciFormatter,
    'sigfig': pticker.SigFigFormatter,
    'simple': pticker.SimpleFormatter,
    'date': mdates.AutoDateFormatter,
    'datestr': mdates.DateFormatter,
    'scalar': mticker.ScalarFormatter,
    'none': mticker.NullFormatter,
    'null': mticker.NullFormatter,
    'func': mticker.FuncFormatter,
    'strmethod': mticker.StrMethodFormatter,
    'formatstr': mticker.FormatStrFormatter,
    'log': mticker.LogFormatterSciNotation,  # NOTE: this is subclass of Mathtext class
    'logit': mticker.LogitFormatter,
    'eng': mticker.EngFormatter,
    'percent': mticker.PercentFormatter,
    'index': pticker._IndexFormatter,
    'e': partial(pticker.FracFormatter, symbol=r'$e$', number=np.e),
    'pi': partial(pticker.FracFormatter, symbol=r'$\pi$', number=np.pi),
    'tau': partial(pticker.FracFormatter, symbol=r'$\tau$', number=2 * np.pi),
    'lat': partial(pticker.SimpleFormatter, negpos='SN'),
    'lon': partial(pticker.SimpleFormatter, negpos='WE', wraprange=(-180, 180)),
    'deg': partial(pticker.SimpleFormatter, suffix='\N{DEGREE SIGN}'),
    'deglat': partial(pticker.SimpleFormatter, negpos='SN', suffix='\N{DEGREE SIGN}'),
    'deglon': partial(pticker.SimpleFormatter, negpos='WE', suffix='\N{DEGREE SIGN}', wraprange=(-180, 180)),  # noqa: E501
    'math': mticker.LogFormatterMathtext,  # deprecated (use SciNotation subclass)
}
if _version_cartopy >= _version('0.18'):
    # NOTE: Will raise error when you try to use these without cartopy >= 0.18
    FORMATTERS['dms'] = partial(pticker._DegreeFormatter, dms=True)
    FORMATTERS['dmslon'] = partial(pticker._LongitudeFormatter, dms=True)
    FORMATTERS['dmslat'] = partial(pticker._LatitudeFormatter, dms=True)
if hasattr(mpolar, 'ThetaFormatter'):
    FORMATTERS['theta'] = mpolar.ThetaFormatter
if hasattr(mdates, 'ConciseDateFormatter'):
    FORMATTERS['concise'] = mdates.ConciseDateFormatter

# The registered scale names and their associated
# `~matplotlib.scale.ScaleBase` classes. See `Scale` for a table.
SCALES = mscale._scale_mapping
SCALE_PRESETS = {
    'quadratic': ('power', 2,),
    'cubic': ('power', 3,),
    'quartic': ('power', 4,),
    'height': ('exp', np.e, -1 / 7, 1013.25, True),
    'pressure': ('exp', np.e, -1 / 7, 1013.25, False),
    'db': ('exp', 10, 1, 0.1, True),
    'idb': ('exp', 10, 1, 0.1, False),
    'np': ('exp', np.e, 1, 1, True),
    'inp': ('exp', np.e, 1, 1, False),
}
mscale.register_scale(pscale.CutoffScale)
mscale.register_scale(pscale.ExpScale)
mscale.register_scale(pscale.LogScale)
mscale.register_scale(pscale.LinearScale)
mscale.register_scale(pscale.LogitScale)
mscale.register_scale(pscale.FuncScale)
mscale.register_scale(pscale.PowerScale)
mscale.register_scale(pscale.SymmetricalLogScale)
mscale.register_scale(pscale.InverseScale)
mscale.register_scale(pscale.SineLatitudeScale)
mscale.register_scale(pscale.MercatorLatitudeScale)

# Default keyword args for `~mpl_toolkits.basemap.Basemap` projections.
# `~mpl_toolkits.basemap` will raise an error if you don't provide them,
# so ProPlot imposes some sensible default behavior.
BASEMAP_KW_DEFAULTS = {
    'eck4': {'lon_0': 0},
    'geos': {'lon_0': 0},
    'hammer': {'lon_0': 0},
    'moll': {'lon_0': 0},
    'kav7': {'lon_0': 0},
    'sinu': {'lon_0': 0},
    'vandg': {'lon_0': 0},
    'mbtfpq': {'lon_0': 0},
    'robin': {'lon_0': 0},
    'ortho': {'lon_0': 0, 'lat_0': 0},
    'nsper': {'lon_0': 0, 'lat_0': 0},
    'aea': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'eqdc': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'cass': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'gnom': {'lon_0': 0, 'lat_0': 90, 'width': 15000e3, 'height': 15000e3},
    'poly': {'lon_0': 0, 'lat_0': 0, 'width': 10000e3, 'height': 10000e3},
    'npaeqd': {'lon_0': 0, 'boundinglat': 10},  # NOTE: everything breaks if you
    'nplaea': {'lon_0': 0, 'boundinglat': 10},  # try to set boundinglat to zero
    'npstere': {'lon_0': 0, 'boundinglat': 10},
    'spaeqd': {'lon_0': 0, 'boundinglat': -10},
    'splaea': {'lon_0': 0, 'boundinglat': -10},
    'spstere': {'lon_0': 0, 'boundinglat': -10},
    'lcc': {
        'lon_0': 0, 'lat_0': 40, 'lat_1': 35, 'lat_2': 45,  # use cartopy defaults
        'width': 20000e3, 'height': 15000e3
    },
    'tmerc': {
        'lon_0': 0, 'lat_0': 0, 'width': 10000e3, 'height': 10000e3
    },
    'merc': {
        'llcrnrlat': -80, 'urcrnrlat': 84, 'llcrnrlon': -180, 'urcrnrlon': 180
    },
    'omerc': {
        'lat_0': 0, 'lon_0': 0, 'lat_1': -10, 'lat_2': 10,
        'lon_1': 0, 'lon_2': 0, 'width': 10000e3, 'height': 10000e3
    },
}
BASEMAP_PROJ_ALIASES = {  # aliases to basemap naming conventions
    'eqc': 'cyl',
    'pcarree': 'cyl',
}

# Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes.
CARTOPY_PROJS = {}
if ccrs is not object:
    CARTOPY_PROJS.update({
        'aitoff': pcrs.Aitoff,
        'hammer': pcrs.Hammer,
        'kav7': pcrs.KavrayskiyVII,
        'wintri': pcrs.WinkelTripel,
        'npgnom': pcrs.NorthPolarGnomonic,
        'spgnom': pcrs.SouthPolarGnomonic,
        'npaeqd': pcrs.NorthPolarAzimuthalEquidistant,
        'spaeqd': pcrs.SouthPolarAzimuthalEquidistant,
        'nplaea': pcrs.NorthPolarLambertAzimuthalEqualArea,
        'splaea': pcrs.SouthPolarLambertAzimuthalEqualArea,
    })
    CARTOPY_PROJS_UNAVAIL = {
        'aea': 'AlbersEqualArea',
        'aeqd': 'AzimuthalEquidistant',
        'cyl': 'PlateCarree',  # only basemap name not matching PROJ
        'eck1': 'EckertI',
        'eck2': 'EckertII',
        'eck3': 'EckertIII',
        'eck4': 'EckertIV',
        'eck5': 'EckertV',
        'eck6': 'EckertVI',
        'eqc': 'PlateCarree',  # actual PROJ name
        'eqdc': 'EquidistantConic',
        'eqearth': 'EqualEarth',  # better looking Robinson; not in basemap
        'euro': 'EuroPP',  # Europe; not in basemap or PROJ
        'geos': 'Geostationary',
        'gnom': 'Gnomonic',
        'igh': 'InterruptedGoodeHomolosine',  # not in basemap
        'laea': 'LambertAzimuthalEqualArea',
        'lcc': 'LambertConformal',
        'lcyl': 'LambertCylindrical',  # not in basemap or PROJ
        'merc': 'Mercator',
        'mill': 'Miller',
        'moll': 'Mollweide',
        'npstere': 'NorthPolarStereo',  # np/sp stuff not in PROJ
        'nsper': 'NearsidePerspective',
        'ortho': 'Orthographic',
        'osgb': 'OSGB',  # UK; not in basemap or PROJ
        'osni': 'OSNI',  # Ireland; not in basemap or PROJ
        'pcarree': 'PlateCarree',  # common alternate name
        'robin': 'Robinson',
        'rotpole': 'RotatedPole',
        'sinu': 'Sinusoidal',
        'spstere': 'SouthPolarStereo',
        'stere': 'Stereographic',
        'tmerc': 'TransverseMercator',
        'utm': 'UTM',  # not in basemap
    }
    for _key, _cls in list(CARTOPY_PROJS_UNAVAIL.items()):
        if hasattr(ccrs, _cls):
            CARTOPY_PROJS[_key] = getattr(ccrs, _cls)
            del CARTOPY_PROJS_UNAVAIL[_key]
    if CARTOPY_PROJS_UNAVAIL:
        warnings._warn_proplot(
            'Cartopy projection(s) '
            + ', '.join(map(repr, CARTOPY_PROJS_UNAVAIL.values()))
            + 'are unavailable. Consider updating to cartopy >= 0.17.0.'
        )
CARTOPY_KW_ALIASES = {  # use PROJ shorthands instead of verbose cartopy names
    'lat_0': 'central_latitude',
    'lon_0': 'central_longitude',
    'lat_min': 'min_latitude',
    'lat_max': 'max_latitude',
}

# Resolution aliases
# NOTE: Maximum basemap resolutions are much finer than cartopy
CARTOPY_RESOS = {
    'lo': '110m',
    'med': '50m',
    'hi': '10m',
    'x-hi': '10m',  # extra high
    'xx-hi': '10m',  # extra extra high
}
BASEMAP_RESOS = {
    'lo': 'c',  # coarse
    'med': 'l',
    'hi': 'i',
    'x-hi': 'h',
    'xx-hi': 'f',  # fine
}

# Geographic feature properties
CARTOPY_FEATURES = {  # positional arguments passed to NaturalEarthFeature
    'land': ('physical', 'land'),
    'ocean': ('physical', 'ocean'),
    'lakes': ('physical', 'lakes'),
    'coast': ('physical', 'coastline'),
    'rivers': ('physical', 'rivers_lake_centerlines'),
    'borders': ('cultural', 'admin_0_boundary_lines_land'),
    'innerborders': ('cultural', 'admin_1_states_provinces_lakes'),
}
BASEMAP_FEATURES = {  # names of relevant basemap methods
    'land': 'fillcontinents',
    'coast': 'drawcoastlines',
    'rivers': 'drawrivers',
    'borders': 'drawcountries',
    'innerborders': 'drawstates',
}


def Colors(*args, **kwargs):
    """
    Pass all arguments to `Cycle` and return the list of colors from
    the resulting `~cycler.Cycler` object.

    See also
    --------
    cycler.Cycler
    Cycle
    """
    cycle = Cycle(*args, **kwargs)
    return list(dict_['color'] for dict_ in cycle)


def _modify_colormap(cmap, *, cut, left, right, reverse, shift, alpha, samples):
    """
    Modify colormap using a variety of methods.
    """
    if cut is not None or left is not None or right is not None:
        if isinstance(cmap, pcolors.ListedColormap):
            if cut is not None:
                warnings._warn_proplot(
                    "Invalid argument 'cut' for ListedColormap. Ignoring."
                )
            cmap = cmap.truncate(left=left, right=right)
        else:
            cmap = cmap.cut(cut, left=left, right=right)
    if reverse:
        cmap = cmap.reversed()
    if shift is not None:
        cmap = cmap.shifted(shift)
    if alpha is not None:
        cmap = cmap.copy(alpha=alpha)
    if samples is not None:
        if isinstance(cmap, pcolors.ListedColormap):
            cmap = cmap.copy(N=samples)
        else:
            cmap = cmap.to_listed(samples)
    return cmap


@warnings._rename_kwargs('v0.7', shade='luminance')
def Colormap(
    *args, name=None, listmode='perceptual', to_listed=False, cycle=None,
    save=False, save_kw=None, **kwargs
):
    """
    Generate, retrieve, modify, and/or merge instances of
    `~proplot.colors.PerceptuallyUniformColormap`,
    `~proplot.colors.LinearSegmentedColormap`, and
    `~proplot.colors.ListedColormap`.
    Used to interpret the `cmap` and `cmap_kw` arguments
    when passed to any plotting method wrapped by
    `~proplot.axes.apply_cmap`.

    Parameters
    ----------
    *args : colormap-spec
        Positional arguments that individually generate colormaps. If more than
        one argument is passed, the resulting colormaps are *merged* with
        `~proplot.colors.LinearSegmentedColormap.append`
        or `~proplot.colors.ListedColormap.append`.
        Arguments are interpreted as follows:

        * If the argument is a `~matplotlib.colors.Colormap` or a registered
          colormap name, nothing more is done.
        * If a filename string with valid extension, the colormap data will
          be loaded. See `~proplot.config.register_cmaps` and
          `~proplot.config.register_cycles`.
        * If RGB tuple or color string, a
          `~proplot.colors.PerceptuallyUniformColormap` is generated with
          `~proplot.colors.PerceptuallyUniformColormap.from_color`. If the
          string ends in ``'_r'``, the monochromatic map will be *reversed*,
          i.e. will go from dark to light instead of light to dark.
        * If list of RGB tuples or color strings, a
          `~proplot.colors.PerceptuallyUniformColormap` is generated with
          `~proplot.colors.PerceptuallyUniformColormap.from_list`.
        * If dictionary, a `~proplot.colors.PerceptuallyUniformColormap` is
          generated with `~proplot.colors.PerceptuallyUniformColormap.from_hsl`.
          The dictionary should contain the keys ``'hue'``, ``'saturation'``,
          ``'luminance'``, and optionally ``'alpha'``, or their aliases (see below).

    name : str, optional
        Name under which the final colormap is registered. It can then be
        reused by passing ``cmap='name'`` to plotting functions like
        `~matplotlib.axes.Axes.contourf`.
    listmode : {'perceptual', 'linear', 'listed'}, optional
        Controls how colormaps are generated when you input list(s) of colors.
        If ``'perceptual'``, a `~proplot.colors.PerceptuallyUniformColormap`
        is generated with `~proplot.colors.PerceptuallyUniformColormap.from_list`.
        If ``'linear'``, a `~matplotlib.colors.LinearSegmentedColormap` is
        generated with `~matplotlib.colors.LinearSegmentedColormap.from_list`.
        If ``'listed'``, a `~matplotlib.colors.ListedColormap` is generated.
        Default is ``'perceptual'`` when calling `Colormap` directly and
        ``'listed'`` when `Colormap` is called by `Cycle`.
    samples : int or list of int, optional
        For `~proplot.colors.LinearSegmentedColormap`\\ s, this is used to
        generate `~proplot.colors.ListedColormap`\\ s with
        `~proplot.colors.LinearSegmentedColormap.to_listed`. For
        `~proplot.colors.ListedColormap`\\ s, this is used to updates the
        number of colors in the cycle. If `samples` is integer, it applies
        to the final *merged* colormap. If it is a list of integers,
        it applies to each input colormap individually.
    to_listed : bool, optional
        If ``True``, when the final colormap is a
        `~proplot.colors.ListedColormap`, we leave it alone, but when it is a
        `~proplot.colors.LinearSegmentedColormap`, we always call
        `~proplot.colors.LinearSegmentedColormap.to_listed` with a
        default `samples` value of ``10``. This argument is not
        necessary if you provide the `samples` argument.
    left, right : float or list of float, optional
        Truncate the left or right edges of the colormap.
        Passed to `~proplot.colors.LinearSegmentedColormap.truncate`.
        If float, these apply to the final *merged* colormap. If list
        of float, these apply to each input colormap individually.
    cut : float or list of float, optional
        Cut out the center of the colormap. Passed to
        `~proplot.colors.LinearSegmentedColormap.cut`. If float,
        this applies to the final *merged* colormap. If list of float,
        these apply to each input colormap individually.
    reverse : bool or list of bool, optional
        Reverse the colormap. Passed to
        `~proplot.colors.LinearSegmentedColormap.reversed`. If
        float, this applies to the final *merged* colormap. If list of
        float, these apply to each input colormap individually.
    shift : float or list of float, optional
        Cyclically shift the colormap.
        Passed to `~proplot.colors.LinearSegmentedColormap.shifted`.
        If float, this applies to the final *merged* colormap. If list of
        float, these apply to each input colormap individually.
    a
        Shorthand for `alpha`.
    alpha, a : channel-spec or list of channel-spec, optional
        The opacity of the colormap or the opacity gradation. Passed to
        `proplot.colors.LinearSegmentedColormap.set_alpha`
        or `proplot.colors.ListedColormap.set_alpha`. If float, this applies
        to the final *merged* colormap. If list of float, these apply to
        each colormap individually.
    h, s, l, c
        Shorthands for `hue`, `luminance`, `saturation`, and `chroma`.
    hue, saturation, luminance : channel-spec or list of channel-spec, optional
        The channel value(s) used to generate colormaps with
        `~proplot.colors.PerceptuallyUniformColormap.from_hsl` and
        `~proplot.colors.PerceptuallyUniformColormap.from_color`.

        * If you provided no positional arguments, these are used to create
          an arbitrary perceptually uniform colormap with
          `~proplot.colors.PerceptuallyUniformColormap.from_hsl`. This
          is an alternative to passing a dictionary as a positional argument
          with `hue`, `saturation`, and `luminance` as dictionary keys (see `args`).
        * If you did provide positional arguments, and any of them are
          color specifications, these control the look of monochromatic colormaps
          generated with `~proplot.colors.PerceptuallyUniformColormap.from_color`.
          To use different values for each colormap, pass a list of floats instead
          of a scalar. Note the default `luminance` is ``90`` if `to_listed`
          is ``True`` and ``100`` otherwise.

    chroma
        Alias for `saturation`.
    cycle : str or list of str, optional
        The registered cycle name or a list of colors used to interpret cycle
        color strings like ``'C0'`` and ``'C2'``. Default is from the active
        property cycler. This lets you make monochromatic colormaps using
        colors selected from arbitrary property cycles.
    save : bool, optional
        Whether to call the colormap/color cycle save method, i.e.
        `proplot.colors.LinearSegmentedColormap.save` or
        `proplot.colors.ListedColormap.save`.
    save_kw : dict-like, optional
        Ignored if `save` is ``False``. Passed to the colormap/color cycle
        save method, i.e. `proplot.colors.LinearSegmentedColormap.save` or
        `proplot.colors.ListedColormap.save`.

    Other parameters
    ----------------
    **kwargs
        Passed to `proplot.colors.LinearSegmentedColormap.copy`,
        `proplot.colors.PerceptuallyUniformColormap.copy`, or
        `proplot.colors.ListedColormap.copy`.

    Returns
    -------
    `~matplotlib.colors.Colormap`
        A `~proplot.colors.LinearSegmentedColormap` or
        `~proplot.colors.ListedColormap` instance.

    See also
    --------
    matplotlib.colors.Colormap
    matplotlib.colors.LinearSegmentedColormap
    matplotlib.colors.ListedColormap
    proplot.axes.apply_cmap
    Cycle
    """
    # Helper function
    # NOTE: Very careful here! Try to support common use cases. For example
    # adding opacity gradations to colormaps with Colormap('cmap', alpha=(0.5, 1))
    # or sampling maps with Colormap('cmap', samples=np.linspace(0, 1, 11)) should
    # be allowable.
    # If *args is singleton try to preserve it.
    def _pop_modification(key):
        value = kwargs.pop(key, None)
        if not np.iterable(value) or isinstance(value, str):
            values = (None,) * len(args)
        elif len(args) == len(value):
            values, value = tuple(values), None
        elif len(args) == 1:  # e.g. Colormap('cmap', alpha=(0.5, 1))
            values = (None,)
        else:
            raise ValueError(
                f'Got {len(args)} colormap-specs '
                f'but {len(value)} values for {key!r}.'
            )
        return value, values

    # Parse keyword args that can apply to the merged colormap or each
    # colormap individually.
    hsla = _pop_props(kwargs, 'hsla')
    if not args and hsla.keys() - {'alpha'}:
        args = (hsla,)
    else:
        kwargs.update(hsla)
    cut, cuts = _pop_modification('cut')
    left, lefts = _pop_modification('left')
    right, rights = _pop_modification('right')
    shift, shifts = _pop_modification('shift')
    reverse, reverses = _pop_modification('reverse')
    samples, sampless = _pop_modification('samples')
    alpha, alphas = _pop_modification('alpha')
    luminance, luminances = _pop_modification('luminance')
    saturation, saturations = _pop_modification('saturation')
    if luminance is not None:
        luminances = (luminance,) * len(args)
    if saturation is not None:
        saturations = (saturation,) * len(args)

    # Parse input args
    # TODO: Play with using "qualitative" colormaps in realistic examples,
    # how to make colormaps cyclic.
    if not args:
        raise ValueError(
            'Colormap() requires either positional arguments '
            "or 'hue', 'chroma', 'saturation', and/or 'luminance' keywords."
        )
    if listmode not in ('listed', 'linear', 'perceptual'):
        raise ValueError(
            f'Invalid listmode={listmode!r}. Options are: '
            "'listed', 'linear', 'perceptual'."
        )

    # Loop through colormaps
    tmp = '_no_name'
    cmaps = []
    for arg, icut, ileft, iright, ireverse, ishift, isamples, iluminance, isaturation, ialpha in zip(  # noqa: E501
        args, cuts, lefts, rights, reverses, shifts, sampless, luminances, saturations, alphas  # noqa: E501
    ):
        # Load registered colormaps and maps on file
        # TODO: Document how 'listmode' also affects loaded files
        if isinstance(arg, str):
            if '.' in arg and os.path.isfile(arg):
                if listmode == 'listed':
                    arg = pcolors.ListedColormap.from_file(arg)
                else:
                    arg = pcolors.LinearSegmentedColormap.from_file(arg)
            else:
                try:
                    arg = pcolors._cmap_database[arg]
                except KeyError:
                    pass

        # Convert matplotlib colormaps to subclasses
        if isinstance(arg, mcolors.Colormap):
            cmap = pcolors._to_proplot_colormap(arg)

        # Dictionary of hue/sat/luminance values or 2-tuples
        elif isinstance(arg, dict):
            cmap = pcolors.PerceptuallyUniformColormap.from_hsl(tmp, **arg)

        # List of color tuples or color strings, i.e. iterable of iterables
        elif (
            not isinstance(arg, str) and np.iterable(arg)
            and all(np.iterable(color) for color in arg)
        ):
            colors = [to_rgba(color, cycle=cycle) for color in arg]
            if listmode == 'listed':
                cmap = pcolors.ListedColormap(colors, tmp)
            elif listmode == 'linear':
                cmap = pcolors.LinearSegmentedColormap.from_list(tmp, colors)
            else:
                cmap = pcolors.PerceptuallyUniformColormap.from_list(tmp, colors)

        # Monochrome colormap from input color
        # NOTE: Do not print color names in error message. Too long to be useful.
        else:
            jreverse = isinstance(arg, str) and arg[-2:] == '_r'
            if jreverse:
                arg = arg[:-2]
            try:
                color = to_rgba(arg, cycle=cycle)
            except (ValueError, TypeError):
                message = f'Invalid colormap, color cycle, or color {arg!r}.'
                if isinstance(arg, str) and arg[:1] != '#':
                    message += (
                        ' Options are: '
                        + ', '.join(map(repr, pcolors._cmap_database)) + '.'
                    )
                raise ValueError(message)
            if to_listed and iluminance is None:
                iluminance = 90
            cmap = pcolors.PerceptuallyUniformColormap.from_color(
                tmp, color, luminance=iluminance, saturation=isaturation
            )
            ireverse = _not_none(ireverse, False)
            ireverse = ireverse ^ jreverse  # xor

        # Modify the colormap
        cmap = _modify_colormap(
            cmap, cut=icut, left=ileft, right=iright,
            reverse=ireverse, shift=ishift, alpha=ialpha, samples=isamples,
        )
        cmaps.append(cmap)

    # Merge the resulting colormaps
    if len(cmaps) > 1:  # more than one map and modify arbitrary properties
        cmap = cmaps[0].append(*cmaps[1:], **kwargs)
    elif kwargs:  # modify arbitrary properties
        cmap = cmaps[0].copy(**kwargs)
    else:
        cmap = cmaps[0]

    # Modify the colormap
    if to_listed and samples is None and isinstance(cmap, pcolors.LinearSegmentedColormap):  # noqa: E501
        samples = 10
    cmap = _modify_colormap(
        cmap, cut=cut, left=left, right=right,
        reverse=reverse, shift=shift, alpha=alpha, samples=samples
    )

    # Initialize
    if not cmap._isinit:
        cmap._init()

    # Register the colormap
    if name is None:
        name = cmap.name  # may have been modified by e.g. .shifted()
    else:
        cmap.name = name
    pcolors._cmap_database[name] = cmap

    # Save the colormap
    if save:
        save_kw = save_kw or {}
        cmap.save(**save_kw)

    return cmap


def Cycle(*args, N=None, samples=None, name=None, **kwargs):
    """
    Generate and merge `~cycler.Cycler` instances in a variety of ways.
    Used to interpret the `cycle` and `cycle_kw` arguments when passed to
    any plotting method wrapped by `~proplot.axes.apply_cycle`.

    If you just want a list of colors instead of a `~cycler.Cycler` instance,
    use the `Colors` function. If you want a `~cycler.Cycler` instance that
    imposes black as the default color and cycles through properties like
    ``linestyle`` instead, call this function without any positional arguments.

    Parameters
    ----------
    *args : colormap-spec or cycle-spec, optional
        Positional arguments control the *colors* in the `~cycler.Cycler`
        object. If more than one argument is passed, the resulting cycles are
        merged. Arguments are interpreted as follows:

        * If a `~cycler.Cycler`, nothing more is done.
        * If a list of RGB tuples or color strings, these colors are used.
        * If a `~proplot.colors.ListedColormap`, colors from the ``colors``
          attribute are used.
        * If a string cycle name, that `~proplot.colors.ListedColormap`
          is looked up and its ``colors`` attribute is used.
        * In all other cases, the argument is passed to `Colormap`, and
          colors from the resulting `~proplot.colors.LinearSegmentedColormap`
          are used. See the `samples` argument.

        If the last positional argument is numeric, it is used for the
        `samples` keyword argument.
    N
        Shorthand for `samples`.
    samples : float or list of float, optional
        For `~proplot.colors.ListedColormap`\\ s, this is the number of
        colors to select. For example, ``Cycle('538', 4)`` returns the first 4
        colors of the ``'538'`` color cycle.

        For `~proplot.colors.LinearSegmentedColormap`\\ s, this is either a
        list of sample coordinates used to draw colors from the map, or an
        integer number of colors to draw. If the latter, the sample coordinates
        are ``np.linspace(0, 1, samples)``. For example, ``Cycle('Reds', 5)``
        divides the ``'Reds'`` colormap into five evenly spaced colors.
    lw, ls, d, a, m, ms, mew, mec, mfc
        Shorthands for the below keywords.
    linewidth, linestyle, dashes, alpha, marker, markersize, markeredgewidth, \
markeredgecolor, markerfacecolor : spec or list of specs, optional
        Lists of `~matplotlib.lines.Line2D` properties that can be added to
        the `~cycler.Cycler` instance. If the lists have unequal length, they
        will be filled to match the length of the longest list.  See
        `~matplotlib.axes.Axes.set_prop_cycle` for more info on cyclers.
        Also see the `line style reference \
<https://matplotlib.org/2.2.5/gallery/lines_bars_and_markers/line_styles_reference.html>`__,
        the `marker reference \
<https://matplotlib.org/stable/gallery/lines_bars_and_markers/marker_reference.html>`__,
        and the `custom dashes reference \
<https://matplotlib.org/stable/gallery/lines_bars_and_markers/line_demo_dash_control.html>`__.
    linewidths, linestyles, dashes, alphas, markers, markersizes, markeredgewidths, \
markeredgecolors, markerfacecolors
        Aliases for the above keywords.

    Other parameters
    ----------------
    **kwargs
        If the input is not already a `~cycler.Cycler` instance, these are passed
        to `Colormap` and used to build the `~proplot.colors.ListedColormap` from
        which the cycler will draw its colors.

    Returns
    -------
    `~cycler.Cycler`
        A cycler instance that can be passed to
        `~matplotlib.axes.Axes.set_prop_cycle`.

    See also
    --------
    cycler.Cycler
    proplot.axes.apply_cycle
    Colormap
    """
    # Parse keyword arguments that rotate through other properties
    # besides color cycles.
    props = _pop_props(kwargs, 'lines')
    nprops = 0
    samples = _not_none(samples=samples, N=N)  # trigger Colormap default
    for key, value in tuple(props.items()):  # permit in-place modification
        if value is None:
            return
        elif not np.iterable(value) or isinstance(value, str):
            value = (value,)
        elif len(value) != len(args):
            nprops = max(nprops, len(value))
        props[key] = list(value)  # ensure mutable list

    # If args is non-empty, means we want color cycle; otherwise is black
    if not args:
        props['color'] = [mcolors.to_rgba('k')]
        if kwargs:
            warnings._warn_proplot(f'Ignoring Cycle() keyword arg(s) {kwargs}.')

    # Merge cycler objects
    elif all(isinstance(arg, cycler.Cycler) for arg in args):
        if kwargs:
            warnings._warn_proplot(f'Ignoring Cycle() keyword arg(s) {kwargs}.')
        if len(args) == 1:
            return args[0]
        else:
            props = {}
            for arg in args:
                for key, value in arg.by_key():
                    if key not in props:
                        props[key] = []
                    props[key].extend(value)
            return cycler.cycler(**props)

    # Get cycler from a colormap
    else:
        # Get the ListedColormap
        if args and isinstance(args[-1], Number):
            args, samples = args[:-1], _not_none(samples_positional=args[-1], samples=samples)  # noqa: #501
        kwargs.setdefault('listmode', 'listed')
        kwargs.setdefault('to_listed', True)  # triggers default 'samples' value
        cmap = Colormap(*args, name=name, samples=samples, **kwargs)
        name = _not_none(name, cmap.name)
        # Add colors to property dict
        nprops = max(nprops, len(cmap.colors))
        props['color'] = [  # save the tupled version!
            tuple(color) if not isinstance(color, str) else color
            for color in cmap.colors
        ]

    # Build cycler, make sure lengths are the same
    for key, value in props.items():
        if len(value) < nprops:  # double back if necessary
            value[:] = [value[i % len(value)] for i in range(nprops)]
    cycle = cycler.cycler(**props)
    cycle.name = _not_none(name, '_no_name')

    return cycle


def Norm(norm, *args, **kwargs):
    """
    Return an arbitrary `~matplotlib.colors.Normalize` instance. Used to
    interpret the `norm` and `norm_kw` arguments when passed to any plotting
    method wrapped by `~proplot.axes.apply_cmap`. See `this tutorial \
<https://matplotlib.org/stable/tutorials/colors/colormapnorms.html>`__
    for more info.

    Parameters
    ----------
    norm : str or `~matplotlib.colors.Normalize`
        The normalizer specification. If a `~matplotlib.colors.Normalize`
        instance already, the input argument is simply returned. Otherwise,
        `norm` should be a string corresponding to one of the "registered"
        colormap normalizers (see below table).

        If `norm` is a list or tuple and the first element is a "registered"
        normalizer name, subsequent elements are passed to the normalizer class
        as positional arguments.

        .. _norm_table:

        ==========================  =====================================
        Key(s)                      Class
        ==========================  =====================================
        ``'null'``, ``'none'``      `~matplotlib.colors.NoNorm`
        ``'diverging'``, ``'div'``  `~proplot.colors.DivergingNorm`
        ``'segmented'``             `~proplot.colors.LinearSegmentedNorm`
        ``'linear'``                `~matplotlib.colors.Normalize`
        ``'log'``                   `~matplotlib.colors.LogNorm`
        ``'power'``                 `~matplotlib.colors.PowerNorm`
        ``'symlog'``                `~matplotlib.colors.SymLogNorm`
        ==========================  =====================================

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.colors.Normalize` initializer.

    Returns
    -------
    `~matplotlib.colors.Normalize`
        A `~matplotlib.colors.Normalize` instance.

    See also
    --------
    matplotlib.colors.Normalize
    proplot.axes.apply_cmap
    proplot.colors.DiscreteNorm
    """
    if isinstance(norm, mcolors.Normalize):
        return norm

    # Pull out extra args
    if np.iterable(norm) and not isinstance(norm, str):
        norm, args = norm[0], (*norm[1:], *args)
    if not isinstance(norm, str):
        raise ValueError(f'Invalid norm name {norm!r}. Must be string.')

    # Get class
    if norm not in NORMS:
        raise ValueError(
            f'Unknown normalizer {norm!r}. Options are: '
            + ', '.join(map(repr, NORMS.keys())) + '.'
        )
    if norm == 'symlog' and not args and 'linthresh' not in kwargs:
        kwargs['linthresh'] = 1  # special case, needs argument
    return NORMS[norm](*args, **kwargs)


def Locator(locator, *args, **kwargs):
    """
    Return a `~matplotlib.ticker.Locator` instance. This function is used to
    interpret the `xlocator`, `xlocator_kw`, `ylocator`, `ylocator_kw`,
    `xminorlocator`, `xminorlocator_kw`, `yminorlocator`, and
    `yminorlocator_kw` arguments when passed to
    `~proplot.axes.CartesianAxes.format`, and the `locator`, `locator_kw`,
    `minorlocator`, and `minorlocator_kw` arguments when passed to colorbar
    methods wrapped by `~proplot.axes.colorbar_extras`.

    Parameters
    ----------
    locator : `~matplotlib.ticker.Locator`, str, float, or list of float
        The axis locator specification, interpreted as follows:

        * If a `~matplotlib.ticker.Locator` instance already, the input
          argument is simply returned.
        * If a list of numbers, these points are ticked. Returns a
          `~matplotlib.ticker.FixedLocator`.
        * If number, this specifies the *step size* between tick locations.
          Returns a `~matplotlib.ticker.MultipleLocator`.

        Otherwise, `locator` should be a string corresponding to one
        of the "registered" locators (see below table). If `locator` is a
        list or tuple and the first element is a "registered" locator name,
        subsequent elements are passed to the locator class as positional
        arguments.

        .. _locator_table:

        =======================  ============================================  =====================================================================================
        Key                      Class                                         Description
        =======================  ============================================  =====================================================================================
        ``'null'``, ``'none'``   `~matplotlib.ticker.NullLocator`              No ticks
        ``'auto'``               `~matplotlib.ticker.AutoLocator`              Major ticks at sensible locations
        ``'minor'``              `~matplotlib.ticker.AutoMinorLocator`         Minor ticks at sensible locations
        ``'date'``               `~matplotlib.dates.AutoDateLocator`           Default tick locations for datetime axes
        ``'fixed'``              `~matplotlib.ticker.FixedLocator`             Ticks at these exact locations
        ``'index'``              `~matplotlib.ticker.IndexLocator`             Ticks on the non-negative integers
        ``'linear'``             `~matplotlib.ticker.LinearLocator`            Exactly ``N`` ticks encompassing axis limits, spaced as ``numpy.linspace(lo, hi, N)``
        ``'log'``                `~matplotlib.ticker.LogLocator`               For log-scale axes
        ``'logminor'``           `~matplotlib.ticker.LogLocator`               For log-scale axes on the 1st through 9th multiples of each power of the base
        ``'logit'``              `~matplotlib.ticker.LogitLocator`             For logit-scale axes
        ``'logitminor'``         `~matplotlib.ticker.LogitLocator`             For logit-scale axes with ``minor=True`` passed to `~matplotlib.ticker.LogitLocator`
        ``'maxn'``               `~matplotlib.ticker.MaxNLocator`              No more than ``N`` ticks at sensible locations
        ``'multiple'``           `~matplotlib.ticker.MultipleLocator`          Ticks every ``N`` step away from zero
        ``'symlog'``             `~matplotlib.ticker.SymmetricalLogLocator`    For symlog-scale axes
        ``'symlogminor'``        `~matplotlib.ticker.SymmetricalLogLocator`    For symlog-scale axes on the 1st through 9th multiples of each power of the base
        ``'theta'``              `~matplotlib.projections.polar.ThetaLocator`  Like the base locator but default locations are every `numpy.pi`/8 radians
        ``'year'``               `~matplotlib.dates.YearLocator`               Ticks every ``N`` years
        ``'month'``              `~matplotlib.dates.MonthLocator`              Ticks every ``N`` months
        ``'weekday'``            `~matplotlib.dates.WeekdayLocator`            Ticks every ``N`` weekdays
        ``'day'``                `~matplotlib.dates.DayLocator`                Ticks every ``N`` days
        ``'hour'``               `~matplotlib.dates.HourLocator`               Ticks every ``N`` hours
        ``'minute'``             `~matplotlib.dates.MinuteLocator`             Ticks every ``N`` minutes
        ``'second'``             `~matplotlib.dates.SecondLocator`             Ticks every ``N`` seconds
        ``'microsecond'``        `~matplotlib.dates.MicrosecondLocator`        Ticks every ``N`` microseconds
        ``'lon'``, ``'deglon'``  `~proplot.ticker.LongitudeLocator`            Longitude gridlines at sensible decimal locations
        ``'lat'``, ``'deglat'``  `~proplot.ticker.LatitudeLocator`             Latitude gridlines at sensible decimal locations
        ``'dms'``                `~proplot.ticker._DegreeLocator`              Gridlines on nice minute and second intervals
        ``'dmslon'``             `~proplot.ticker.LongitudeLocator`            Longitude gridlines on nice minute and second intervals
        ``'dmslat'``             `~proplot.ticker.LatitudeLocator`             Latitude gridlines on nice minute and second intervals
        =======================  ============================================  =====================================================================================

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Locator` class.

    Returns
    -------
    `~matplotlib.ticker.Locator`
        A `~matplotlib.ticker.Locator` instance.

    See also
    --------
    matplotlib.ticker.Locator
    proplot.axes.CartesianAxes.format
    proplot.axes.PolarAxes.format
    proplot.axes.GeoAxes.format
    proplot.axes.colorbar_extras
    Formatter
    """  # noqa
    if isinstance(locator, mticker.Locator):
        return locator

    # Pull out extra args
    if np.iterable(locator) and not isinstance(locator, str) and not all(
        isinstance(num, Number) for num in locator
    ):
        locator, args = locator[0], (*locator[1:], *args)

    # Get the locator
    if isinstance(locator, str):  # dictionary lookup
        # Shorthands and defaults
        if locator in ('logminor', 'logitminor', 'symlogminor'):
            locator, _ = locator.split('minor')
            if locator == 'logit':
                kwargs.setdefault('minor', True)
            else:
                kwargs.setdefault('subs', np.arange(1, 10))
        elif locator == 'index':
            args = args or (1,)
            if len(args) == 1:
                args = (*args, 0)
        # Lookup
        if locator not in LOCATORS:
            raise ValueError(
                f'Unknown locator {locator!r}. Options are '
                + ', '.join(map(repr, LOCATORS.keys())) + '.'
            )
        locator = LOCATORS[locator](*args, **kwargs)
    elif isinstance(locator, Number):  # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    elif np.iterable(locator):
        locator = mticker.FixedLocator(np.sort(locator), *args, **kwargs)
    else:
        raise ValueError(f'Invalid locator {locator!r}.')
    return locator


def Formatter(formatter, *args, date=False, index=False, **kwargs):
    """
    Return a `~matplotlib.ticker.Formatter` instance. This function is used to
    interpret the `xformatter`, `xformatter_kw`, `yformatter`, and
    `yformatter_kw` arguments when passed to
    `~proplot.axes.CartesianAxes.format`, and the `formatter`
    and `formatter_kw` arguments when passed to colorbar methods wrapped by
    `~proplot.axes.colorbar_extras`.

    Parameters
    ----------
    formatter : `~matplotlib.ticker.Formatter`, str, list of str, or function
        The axis formatter specification, interpreted as follows:

        * If a `~matplotlib.ticker.Formatter` instance already, the input
          argument is simply returned.
        * If list of strings, the ticks are labeled with these strings. Returns
          a `~matplotlib.ticker.FixedFormatter` if `index` is ``False`` or an
          `~matplotlib.ticker.IndexFormatter` if `index` is ``True``.
        * If a function, the labels will be generated using this function.
          Returns a `~matplotlib.ticker.FuncFormatter`.
        * If a string containing ``{x}`` or ``{x:...}``, ticks will be
          formatted by calling ``string.format(x=number)``.
        * If a string containing ``'%'`` and `date` is ``False``, ticks will be
          formatted using the C-style ``string % number`` method. See
          `this page <https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting>`__
          for a review. Returns a `~matplotlib.ticker.FormatStrFormatter`.
        * If a string containing ``'%'`` and `date` is ``True``, *datetime*
          `string % number`` formatting is used. See
          `this page <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`__
          for a review. Returns a `~matplotlib.ticker.DateFormatter`.

        Otherwise, `formatter` should be a string corresponding to one of the
        "registered" formatters or formatter presets (see below table). If
        `formatter` is a list or tuple and the first element is a "registered"
        formatter name, subsequent elements are passed to the formatter class
        as positional arguments.

        .. _tau: https://tauday.com/tau-manifesto

        .. _formatter_table:

        ======================  ==============================================  =================================================================
        Key                     Class                                           Description
        ======================  ==============================================  =================================================================
        ``'null'``, ``'none'``  `~matplotlib.ticker.NullFormatter`              No tick labels
        ``'auto'``              `~proplot.ticker.AutoFormatter`                 New default tick labels for axes
        ``'sci'``               `~proplot.ticker.SciFormatter`                  Format ticks with scientific notation.
        ``'simple'``            `~proplot.ticker.SimpleFormatter`               New default tick labels for e.g. contour labels
        ``'sigfig'``            `~proplot.ticker.SigFigFormatter`               Format labels using the first ``N`` significant digits
        ``'frac'``              `~proplot.ticker.FracFormatter`                 Rational fractions
        ``'date'``              `~matplotlib.dates.AutoDateFormatter`           Default tick labels for datetime axes
        ``'concise'``           `~matplotlib.dates.ConciseDateFormatter`        More concise date labels introduced in matplotlib 3.1
        ``'datestr'``           `~matplotlib.dates.DateFormatter`               Date formatting with C-style ``string % format`` notation
        ``'eng'``               `~matplotlib.ticker.EngFormatter`               Engineering notation
        ``'fixed'``             `~matplotlib.ticker.FixedFormatter`             List of strings
        ``'formatstr'``         `~matplotlib.ticker.FormatStrFormatter`         From C-style ``string % format`` notation
        ``'func'``              `~matplotlib.ticker.FuncFormatter`              Use an arbitrary function
        ``'index'``             `~matplotlib.ticker.IndexFormatter`             List of strings corresponding to non-negative integer positions
        ``'log'``               `~matplotlib.ticker.LogFormatterSciNotation`    For log-scale axes with scientific notation
        ``'logit'``             `~matplotlib.ticker.LogitFormatter`             For logistic-scale axes
        ``'percent'``           `~matplotlib.ticker.PercentFormatter`           Trailing percent sign
        ``'scalar'``            `~matplotlib.ticker.ScalarFormatter`            Old default tick labels for axes
        ``'strmethod'``         `~matplotlib.ticker.StrMethodFormatter`         From the ``string.format`` method
        ``'theta'``             `~matplotlib.projections.polar.ThetaFormatter`  Formats radians as degrees, with a degree symbol
        ``'e'``                 `~proplot.ticker.FracFormatter` preset          Fractions of *e*
        ``'pi'``                `~proplot.ticker.FracFormatter` preset          Fractions of :math:`\\pi`
        ``'tau'``               `~proplot.ticker.FracFormatter` preset          Fractions of the `one true circle constant <tau_>`_ :math:`\\tau`
        ``'lat'``               `~proplot.ticker.AutoFormatter` preset          Cardinal "SN" indicator
        ``'lon'``               `~proplot.ticker.AutoFormatter` preset          Cardinal "WE" indicator
        ``'deg'``               `~proplot.ticker.AutoFormatter` preset          Trailing degree symbol
        ``'deglat'``            `~proplot.ticker.AutoFormatter` preset          Trailing degree symbol and cardinal "SN" indicator
        ``'deglon'``            `~proplot.ticker.AutoFormatter` preset          Trailing degree symbol and cardinal "WE" indicator
        ``'dms'``               `~cartopy.mpl.ticker._PlateCarreeFormatter`     Labels with degree/minute/second support
        ``'dmslon'``            `~cartopy.mpl.ticker.LongitudeFormatter`        Longitude labels with degree/minute/second support
        ``'dmslat'``            `~cartopy.mpl.ticker.LatitudeFormatter`         Latitude labels with degree/minute/second support
        ======================  ==============================================  =================================================================

    date : bool, optional
        Toggles the behavior when `formatter` contains a ``'%'`` sign (see
        above).
    index : bool, optional
        Controls the behavior when `formatter` is a list of strings (see
        above).

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Formatter` class.

    Returns
    -------
    `~matplotlib.ticker.Formatter`
        A `~matplotlib.ticker.Formatter` instance.

    See also
    --------
    matplotlib.ticker.Formatter
    proplot.axes.CartesianAxes.format
    proplot.axes.PolarAxes.format
    proplot.axes.GeoAxes.format
    proplot.axes.colorbar_extras
    Locator
    """  # noqa
    if isinstance(formatter, mticker.Formatter):  # formatter object
        return formatter

    # Pull out extra args
    if np.iterable(formatter) and not isinstance(formatter, str) and not all(
        isinstance(item, str) for item in formatter
    ):
        formatter, args = formatter[0], (*formatter[1:], *args)

    # Get the formatter
    if isinstance(formatter, str):  # assumption is list of strings
        # Format strings
        if re.search(r'{x(:.+)?}', formatter):
            # string.format() formatting
            formatter = mticker.StrMethodFormatter(
                formatter, *args, **kwargs
            )
        elif '%' in formatter:
            # %-style formatting
            if date:
                formatter = mdates.DateFormatter(
                    formatter, *args, **kwargs
                )
            else:
                formatter = mticker.FormatStrFormatter(
                    formatter, *args, **kwargs
                )
        elif formatter in FORMATTERS:
            # Lookup
            formatter = FORMATTERS[formatter](*args, **kwargs)
        else:
            raise ValueError(
                f'Unknown formatter {formatter!r}. Options are '
                + ', '.join(map(repr, FORMATTERS.keys())) + '.'
            )
    elif callable(formatter):
        # Function
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif np.iterable(formatter):
        # List of strings
        if index:
            formatter = pticker._IndexFormatter(formatter)
        else:
            formatter = mticker.FixedFormatter(formatter)
    else:
        raise ValueError(f'Invalid formatter {formatter!r}.')
    return formatter


def Scale(scale, *args, **kwargs):
    """
    Return a `~matplotlib.scale.ScaleBase` instance. This function is used to
    interpret the `xscale`, `xscale_kw`, `yscale`, and `yscale_kw` arguments
    when passed to `~proplot.axes.CartesianAxes.format`.

    Parameters
    ----------
    scale : `~matplotlib.scale.ScaleBase`, str, or (str, ...)
        The axis scale specification. If a `~matplotlib.scale.ScaleBase`
        instance already, the input argument is simply returned. Otherwise,
        `scale` should be a string corresponding to one of the
        "registered" axis scales or axis scale presets (see below table).

        If `scale` is a list or tuple and the first element is a
        "registered" scale name, subsequent elements are passed to the
        scale class as positional arguments.

        .. _scale_table:

        =================  ======================================  ===============================================
        Key                Class                                   Description
        =================  ======================================  ===============================================
        ``'linear'``       `~proplot.scale.LinearScale`            Linear
        ``'log'``          `~proplot.scale.LogScale`               Logarithmic
        ``'symlog'``       `~proplot.scale.SymmetricalLogScale`    Logarithmic beyond finite space around zero
        ``'logit'``        `~proplot.scale.LogitScale`             Logistic
        ``'inverse'``      `~proplot.scale.InverseScale`           Inverse
        ``'function'``     `~proplot.scale.FuncScale`              Arbitrary forward and backwards transformations
        ``'sine'``         `~proplot.scale.SineLatitudeScale`      Sine function (in degrees)
        ``'mercator'``     `~proplot.scale.MercatorLatitudeScale`  Mercator latitude function (in degrees)
        ``'exp'``          `~proplot.scale.ExpScale`               Arbitrary exponential function
        ``'power'``        `~proplot.scale.PowerScale`             Arbitrary power function
        ``'cutoff'``       `~proplot.scale.CutoffScale`            Arbitrary piecewise linear transformations
        ``'quadratic'``    `~proplot.scale.PowerScale` (preset)    Quadratic function
        ``'cubic'``        `~proplot.scale.PowerScale` (preset)    Cubic function
        ``'quartic'``      `~proplot.scale.PowerScale` (preset)    Quartic function
        ``'db'``           `~proplot.scale.ExpScale` (preset)      Ratio expressed as `decibels <db_>`_
        ``'np'``           `~proplot.scale.ExpScale` (preset)      Ratio expressed as `nepers <np_>`_
        ``'idb'``          `~proplot.scale.ExpScale` (preset)      `Decibels <db_>`_ expressed as ratio
        ``'inp'``          `~proplot.scale.ExpScale` (preset)      `Nepers <np_>`_ expressed as ratio
        ``'pressure'``     `~proplot.scale.ExpScale` (preset)      Height (in km) expressed linear in pressure
        ``'height'``       `~proplot.scale.ExpScale` (preset)      Pressure (in hPa) expressed linear in height
        =================  ======================================  ===============================================

        .. _db: https://en.wikipedia.org/wiki/Decibel
        .. _np: https://en.wikipedia.org/wiki/Neper

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.scale.ScaleBase` class.

    Returns
    -------
    `~matplotlib.scale.ScaleBase`
        The scale instance.

    See also
    --------
    matplotlib.scale.ScaleBase
    proplot.axes.CartesianAxes.format
    proplot.axes.CartesianAxes.dualx
    proplot.axes.CartesianAxes.dualy
    """  # noqa
    # NOTE: Why not try to interpret FuncScale arguments, like when lists
    # of numbers are passed to Locator? Because FuncScale *itself* accepts
    # ScaleBase classes as arguments... but constructor functions cannot
    # do anything but return the class instance upon receiving one.
    if isinstance(scale, mscale.ScaleBase):
        return scale

    # Pull out extra args
    if np.iterable(scale) and not isinstance(scale, str):
        scale, args = scale[0], (*scale[1:], *args)
    if not isinstance(scale, str):
        raise ValueError(f'Invalid scale name {scale!r}. Must be string.')

    # Get scale preset
    if scale in SCALE_PRESETS:
        if args or kwargs:
            warnings._warn_proplot(
                f'Scale {scale!r} is a scale *preset*. Ignoring positional '
                'argument(s): {args} and keyword argument(s): {kwargs}. '
            )
        scale, *args = SCALE_PRESETS[scale]

    # Get scale
    scale = scale.lower()
    if scale in SCALES:
        scale = SCALES[scale]
    else:
        raise ValueError(
            f'Unknown scale or preset {scale!r}. Options are '
            + ', '.join(map(repr, list(SCALES) + list(SCALE_PRESETS))) + '.'
        )
    return scale(*args, **kwargs)


def Proj(name, basemap=None, **kwargs):
    """
    Return a `cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap`
    instance. Used to interpret the `proj` and `proj_kw` arguments when
    passed to `~proplot.ui.subplots`.

    Parameters
    ----------
    name : str, `cartopy.crs.Projection`, or `~mpl_toolkits.basemap.Basemap`
        The projection name or projection class instance. If the latter, it
        is simply returned. If the former, it must correspond to one of the
        `PROJ <https://proj.org>`__ projection name shorthands, like in
        basemap.

        The following table lists the valid projection name shorthands, their
        full names (with links to the relevant
        `PROJ documentation <https://proj4.org/operations/projections>`__),
        and whether they are available in the cartopy and basemap packages.
        (added) indicates a projection class that ProPlot has "added"
        to cartopy using the cartopy API.

        .. _proj_table:

        =============  ===============================================  =========  =======
        Key            Name                                             Cartopy    Basemap
        =============  ===============================================  =========  =======
        ``'aea'``      `Albers Equal Area <aea_>`_                      ✓          ✓
        ``'aeqd'``     `Azimuthal Equidistant <aeqd_>`_                 ✓          ✓
        ``'aitoff'``   `Aitoff <aitoff_>`_                              ✓ (added)  ✗
        ``'cass'``     `Cassini-Soldner <cass_>`_                       ✗          ✓
        ``'cea'``      `Cylindrical Equal Area <cea_>`_                 ✗          ✓
        ``'cyl'``      `Cylindrical Equidistant <eqc_>`_                ✓          ✓
        ``'eck1'``     `Eckert I <eck1_>`_                              ✓          ✗
        ``'eck2'``     `Eckert II <eck2_>`_                             ✓          ✗
        ``'eck3'``     `Eckert III <eck3_>`_                            ✓          ✗
        ``'eck4'``     `Eckert IV <eck4_>`_                             ✓          ✓
        ``'eck5'``     `Eckert V <eck5_>`_                              ✓          ✗
        ``'eck6'``     `Eckert VI <eck6_>`_                             ✓          ✗
        ``'eqdc'``     `Equidistant Conic <eqdc_>`_                     ✓          ✓
        ``'eqc'``      `Cylindrical Equidistant <eqc_>`_                ✓          ✓
        ``'eqearth'``  `Equal Earth <eqearth_>`_                        ✓          ✗
        ``'europp'``   Euro PP (Europe)                                 ✓          ✗
        ``'gall'``     `Gall Stereographic Cylindrical <gall_>`_        ✗          ✓
        ``'geos'``     `Geostationary <geos_>`_                         ✓          ✓
        ``'gnom'``     `Gnomonic <gnom_>`_                              ✓          ✓
        ``'hammer'``   `Hammer <hammer_>`_                              ✓ (added)  ✓
        ``'igh'``      `Interrupted Goode Homolosine <igh_>`_           ✓          ✗
        ``'kav7'``     `Kavrayskiy VII <kav7_>`_                        ✓ (added)  ✓
        ``'laea'``     `Lambert Azimuthal Equal Area <laea_>`_          ✓          ✓
        ``'lcc'``      `Lambert Conformal <lcc_>`_                      ✓          ✓
        ``'lcyl'``     Lambert Cylindrical                              ✓          ✗
        ``'mbtfpq'``   `McBryde-Thomas Flat-Polar Quartic <mbtfpq_>`_   ✗          ✓
        ``'merc'``     `Mercator <merc_>`_                              ✓          ✓
        ``'mill'``     `Miller Cylindrical <mill_>`_                    ✓          ✓
        ``'moll'``     `Mollweide <moll_>`_                             ✓          ✓
        ``'npaeqd'``   North-Polar Azimuthal Equidistant                ✓ (added)  ✓
        ``'npgnom'``   North-Polar Gnomonic                             ✓ (added)  ✗
        ``'nplaea'``   North-Polar Lambert Azimuthal                    ✓ (added)  ✓
        ``'npstere'``  North-Polar Stereographic                        ✓          ✓
        ``'nsper'``    `Near-Sided Perspective <nsper_>`_               ✓          ✓
        ``'osni'``     OSNI (Ireland)                                   ✓          ✗
        ``'osgb'``     OSGB (UK)                                        ✓          ✗
        ``'omerc'``    `Oblique Mercator <omerc_>`_                     ✗          ✓
        ``'ortho'``    `Orthographic <ortho_>`_                         ✓          ✓
        ``'pcarree'``  `Cylindrical Equidistant <eqc_>`_                ✓          ✓
        ``'poly'``     `Polyconic <poly_>`_                             ✗          ✓
        ``'rotpole'``  Rotated Pole                                     ✓          ✓
        ``'sinu'``     `Sinusoidal <sinu_>`_                            ✓          ✓
        ``'spaeqd'``   South-Polar Azimuthal Equidistant                ✓ (added)  ✓
        ``'spgnom'``   South-Polar Gnomonic                             ✓ (added)  ✗
        ``'splaea'``   South-Polar Lambert Azimuthal                    ✓ (added)  ✓
        ``'spstere'``  South-Polar Stereographic                        ✓          ✓
        ``'stere'``    `Stereographic <stere_>`_                        ✓          ✓
        ``'tmerc'``    `Transverse Mercator <tmerc_>`_                  ✓          ✓
        ``'utm'``      `Universal Transverse Mercator <utm_>`_          ✓          ✗
        ``'vandg'``    `van der Grinten <vandg_>`_                      ✗          ✓
        ``'wintri'``   `Winkel tripel <wintri_>`_                       ✓ (added)  ✗
        =============  ===============================================  =========  =======

    basemap : bool, optional
        Whether to use the basemap package as opposed to the cartopy package.
        Default is :rc:`basemap`.
    lonlim : 2-tuple of float, optional
        Alternative way to specify `llcrnrlon` and `urcrnrlon` for basemap
        projections.
    latlim : 2-tuple of float, optional
        Alternative way to specify `llcrnrlat` and `urcrnrlat` for basemap
        projections.

    Other parameters
    ----------------
    **kwargs
        Passed to the `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class. For cartopy axes,
        ProPlot translates `lon_0` and `lat_0` to `central_longitude` and
        `central_latitude`.

    Returns
    -------
    proj : `~mpl_toolkits.basemap.Basemap` or `~cartopy.crs.Projection`
        The projection instance.

    See also
    --------
    mpl_toolkits.basemap.Basemap
    cartopy.crs.Projection
    proplot.ui.subplots
    proplot.axes.GeoAxes
    proplot.axes.CartopyAxes
    proplot.axes.BasemapAxes

    References
    ----------
    For more information on map projections, see the
    `wikipedia page <https://en.wikipedia.org/wiki/Map_projection>`__ and the
    `PROJ <https://proj.org>`__ documentation.

    .. _aea: https://proj4.org/operations/projections/aea.html
    .. _aeqd: https://proj4.org/operations/projections/aeqd.html
    .. _aitoff: https://proj4.org/operations/projections/aitoff.html
    .. _cass: https://proj4.org/operations/projections/cass.html
    .. _cea: https://proj4.org/operations/projections/cea.html
    .. _eqc: https://proj4.org/operations/projections/eqc.html
    .. _eck1: https://proj4.org/operations/projections/eck1.html
    .. _eck2: https://proj4.org/operations/projections/eck2.html
    .. _eck3: https://proj4.org/operations/projections/eck3.html
    .. _eck4: https://proj4.org/operations/projections/eck4.html
    .. _eck5: https://proj4.org/operations/projections/eck5.html
    .. _eck6: https://proj4.org/operations/projections/eck6.html
    .. _eqdc: https://proj4.org/operations/projections/eqdc.html
    .. _eqc: https://proj4.org/operations/projections/eqc.html
    .. _eqearth: https://proj4.org/operations/projections/eqearth.html
    .. _gall: https://proj4.org/operations/projections/gall.html
    .. _geos: https://proj4.org/operations/projections/geos.html
    .. _gnom: https://proj4.org/operations/projections/gnom.html
    .. _hammer: https://proj4.org/operations/projections/hammer.html
    .. _igh: https://proj4.org/operations/projections/igh.html
    .. _kav7: https://proj4.org/operations/projections/kav7.html
    .. _laea: https://proj4.org/operations/projections/laea.html
    .. _lcc: https://proj4.org/operations/projections/lcc.html
    .. _mbtfpq: https://proj4.org/operations/projections/mbtfpq.html
    .. _merc: https://proj4.org/operations/projections/merc.html
    .. _mill: https://proj4.org/operations/projections/mill.html
    .. _moll: https://proj4.org/operations/projections/moll.html
    .. _nsper: https://proj4.org/operations/projections/nsper.html
    .. _omerc: https://proj4.org/operations/projections/omerc.html
    .. _ortho: https://proj4.org/operations/projections/ortho.html
    .. _eqc: https://proj4.org/operations/projections/eqc.html
    .. _poly: https://proj4.org/operations/projections/poly.html
    .. _sinu: https://proj4.org/operations/projections/sinu.html
    .. _stere: https://proj4.org/operations/projections/stere.html
    .. _tmerc: https://proj4.org/operations/projections/tmerc.html
    .. _utm: https://proj4.org/operations/projections/utm.html
    .. _vandg: https://proj4.org/operations/projections/vandg.html
    .. _wintri: https://proj4.org/operations/projections/wintri.html
    """  # noqa
    # Class instances
    is_crs = CRS is not object and isinstance(name, CRS)
    is_basemap = Basemap is not object and isinstance(name, Basemap)
    if is_crs or is_basemap:
        proj = name
        proj._proj_package = 'cartopy' if is_crs else 'basemap'
        if basemap is not None:
            kwargs['basemap'] = basemap
        if kwargs:
            warnings._warn_proplot(f'Ignoring Proj() keyword arg(s): {kwargs!r}.')

    # Invalid
    elif not isinstance(name, str):
        raise ValueError(
            f'Unexpected Proj() argument {name!r}. '
            'Must be name, mpl_toolkits.basemap.Basemap instance, '
            'or cartopy.crs.CRS instance.'
        )

    # Basemap
    elif basemap or basemap is None and rc['basemap']:
        # NOTE: Known issue that basemap sometimes produces backwards maps:
        # https://stackoverflow.com/q/56299971/4970632
        # NOTE: We set rsphere to fix non-conda installed basemap issue:
        # https://github.com/matplotlib/basemap/issues/361
        # NOTE: Unlike cartopy, basemap resolution is configured on
        # initialization and controls *all* features.
        import mpl_toolkits.basemap as mbasemap
        if _version_mpl >= _version('3.3'):
            raise RuntimeError(
                'Basemap is no longer maintained and is incompatible with '
                'matplotlib >= 3.3. Please use cartopy as your cartographic '
                'plotting backend or downgrade to matplotlib <= 3.2.'
            )
        if 'lonlim' in kwargs:
            kwargs['llcrnrlon'], kwargs['urcrnrlon'] = kwargs.pop('lonlim')
        if 'latlim' in kwargs:
            kwargs['llcrnrlat'], kwargs['urcrnrlat'] = kwargs.pop('latlim')
        name = BASEMAP_PROJ_ALIASES.get(name, name)
        kwproj = BASEMAP_KW_DEFAULTS.get(name, {}).copy()
        kwproj.update(kwargs)
        kwproj.setdefault('fix_aspect', True)
        if kwproj.get('lon_0', 0) > 0:
            # Fix issues with Robinson (and related?) projections
            # See: https://stackoverflow.com/questions/56299971/
            # Get both this issue *and* 'no room for axes' issue
            kwproj['lon_0'] -= 360
        if name[:2] in ('np', 'sp'):
            kwproj.setdefault('round', True)
        if name == 'geos':
            kwproj.setdefault('rsphere', (6378137.00, 6356752.3142))
        reso = _not_none(
            reso=kwproj.pop('reso', None),
            resolution=kwproj.pop('resolution', None),
            default=rc['reso']
        )
        try:
            reso = BASEMAP_RESOS[reso]
        except KeyError:
            raise ValueError(
                f'Invalid resolution {reso!r}. Options are: '
                + ', '.join(map(repr, BASEMAP_RESOS)) + '.'
            )
        kwproj.update({'resolution': reso, 'projection': name})
        proj = mbasemap.Basemap(**kwproj)
        proj._proj_package = 'basemap'

    # Cartopy
    else:
        import cartopy.crs  # noqa: F401
        kwproj = {
            CARTOPY_KW_ALIASES.get(key, key): value
            for key, value in kwargs.items()
        }
        crs = CARTOPY_PROJS.get(name, None)
        if name == 'geos':  # fix common mistake
            kwproj.pop('central_latitude', None)
        if 'boundinglat' in kwproj:
            raise ValueError(
                '"boundinglat" must be passed to the ax.format() command '
                'for cartopy axes.'
            )
        if crs is None:
            raise ValueError(
                f'Unknown projection {name!r}. Options are: '
                + ', '.join(map(repr, CARTOPY_PROJS.keys())) + '.'
            )
        proj = crs(**kwproj)
        proj._proj_package = 'cartopy'

    return proj
