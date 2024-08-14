#!/usr/bin/env python3
"""
Utilities for processing input data passed to plotting commands.
"""
import functools
import sys

import numpy as np
import numpy.ma as ma

from . import ic  # noqa: F401
from . import _not_none, warnings

try:
    from cartopy.crs import PlateCarree
except ModuleNotFoundError:
    PlateCarree = object


# Constants
BASEMAP_FUNCS = (  # default latlon=True
    "barbs",
    "contour",
    "contourf",
    "hexbin",
    "imshow",
    "pcolor",
    "pcolormesh",
    "plot",
    "quiver",
    "scatter",
    "streamplot",
    "step",
)
CARTOPY_FUNCS = (  # default transform=PlateCarree()
    "barbs",
    "contour",
    "contourf",
    "fill",
    "fill_between",
    "fill_betweenx",  # NOTE: not sure if these work
    "imshow",
    "pcolor",
    "pcolormesh",
    "plot",
    "quiver",
    "scatter",
    "streamplot",
    "step",
    "tricontour",
    "tricontourf",
    "tripcolor",  # NOTE: not sure why these work
)


def _load_objects():
    """
    Load array-like objects.
    """
    # NOTE: We just want to detect if *input arrays* belong to these types -- and if
    # this is the case, it means the module has already been imported! So, we only
    # try loading these classes within autoformat calls. This saves >500ms of import
    # time. We use ndarray as the default value for unimported types and in loops we
    # are careful to check membership to np.ndarray before anything else.
    global ndarray, DataArray, DataFrame, Series, Index, Quantity
    ndarray = np.ndarray
    DataArray = getattr(sys.modules.get("xarray", None), "DataArray", ndarray)
    DataFrame = getattr(sys.modules.get("pandas", None), "DataFrame", ndarray)
    Series = getattr(sys.modules.get("pandas", None), "Series", ndarray)
    Index = getattr(sys.modules.get("pandas", None), "Index", ndarray)
    Quantity = getattr(sys.modules.get("pint", None), "Quantity", ndarray)


_load_objects()


# Type utilities
def _is_numeric(data):
    """
    Test whether input is numeric array rather than datetime or strings.
    """
    array = _to_numpy_array(data)
    return len(data) and (
        np.issubdtype(array.dtype, np.number)
        or np.issubdtype(array.dtype, object)
        and all(isinstance(_, np.number) for _ in array.flat)
    )


def _is_categorical(data):
    """
    Test whether input is array of strings.
    """
    array = _to_numpy_array(data)
    return len(data) and (
        np.issubdtype(array.dtype, str)
        or np.issubdtype(array.dtype, object)
        and any(isinstance(_, str) for _ in array.flat)
    )


def _is_descending(data):
    """
    Test whether the input data is descending. This is used for auto axis reversal.
    """
    # NOTE: Want this to work with e.g. datetime object arrays and numpy datetime
    # arrays so use try except clause.
    data = _to_numpy_array(data)
    if data.ndim > 1 or data.size < 2:
        return False
    try:
        return all(x != abs(x) for x in np.diff(data))
    except TypeError:
        return False


def _to_duck_array(data, strip_units=False):
    """
    Convert arbitrary input to duck array. Preserve array containers with metadata.
    """
    _load_objects()
    if data is None:
        raise ValueError("Invalid data None.")
    if not isinstance(
        data, (ndarray, DataArray, DataFrame, Series, Index, Quantity)
    ) or not np.iterable(data):
        # WARNING: this strips e.g. scalar DataArray metadata
        data = _to_numpy_array(data)
    if strip_units:  # used for z coordinates that cannot have units
        if isinstance(data, (ndarray, Quantity)):
            if Quantity is not ndarray and isinstance(data, Quantity):
                data = data.magnitude
        elif isinstance(data, DataArray):
            if Quantity is not ndarray and isinstance(data.data, Quantity):
                data = data.copy(deep=False)
                data.data = data.data.magnitude
    return data


def _to_numpy_array(data, strip_units=False):
    """
    Convert arbitrary input to numpy array. Preserve masked arrays and unit arrays.
    """
    _load_objects()
    if data is None:
        raise ValueError("Invalid data None.")
    if isinstance(data, ndarray):
        pass
    elif isinstance(data, DataArray):
        data = data.data  # support pint quantities that get unit-stripped later
    elif isinstance(data, (DataFrame, Series, Index)):
        data = data.values
    if Quantity is not ndarray and isinstance(data, Quantity):
        if strip_units:
            return np.atleast_1d(data.magnitude)
        else:
            return np.atleast_1d(data.magnitude) * data.units
    else:
        return np.atleast_1d(data)  # natively preserves masked arrays


def _to_masked_array(data, *, copy=False):
    """
    Convert numpy array to masked array with consideration for datetimes and quantities.
    """
    units = None
    if ndarray is not Quantity and isinstance(data, Quantity):
        data, units = data.magnitude, data.units
    else:
        data = _to_numpy_array(data)
    if data.dtype == "O":
        data = ma.array(data, mask=False)
    else:
        data = ma.masked_invalid(data, copy=copy)
    if np.issubdtype(data.dtype, np.integer):
        data = data.astype(np.float64)
    if np.issubdtype(data.dtype, np.number):
        data.fill_value *= np.nan  # default float fill_value is 1e+20 or 1e+20 + 0j
    else:
        pass  # leave with default fill_value (e.g. NaT for datetime data)
    return data, units


# Input data transformations
def _to_edges(x, y, z):
    """
    Enforce that coordinates are edges. Convert from centers if possible.
    """
    from ..utils import edges, edges2d

    xlen, ylen = x.shape[-1], y.shape[0]
    if z.ndim == 2 and z.shape[1] == xlen and z.shape[0] == ylen:
        # Get edges given centers
        if all(z.ndim == 1 and z.size > 1 and _is_numeric(z) for z in (x, y)):
            x = edges(x)
            y = edges(y)
        else:
            if x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1 and _is_numeric(x):
                x = edges2d(x)
            if y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1 and _is_numeric(y):
                y = edges2d(y)
    elif z.shape[-1] != xlen - 1 or z.shape[0] != ylen - 1:
        # Helpful error message
        raise ValueError(
            f"Input shapes x {x.shape} and y {y.shape} must match "
            f"array centers {z.shape} or "
            f"array borders {tuple(i + 1 for i in z.shape)}."
        )
    return x, y


def _to_centers(x, y, z):
    """
    Enforce that coordinates are centers. Convert from edges if possible.
    """
    xlen, ylen = x.shape[-1], y.shape[0]
    if z.ndim == 2 and z.shape[1] == xlen - 1 and z.shape[0] == ylen - 1:
        # Get centers given edges
        if all(z.ndim == 1 and z.size > 1 and _is_numeric(z) for z in (x, y)):
            x = 0.5 * (x[1:] + x[:-1])
            y = 0.5 * (y[1:] + y[:-1])
        else:
            if x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1 and _is_numeric(x):
                x = 0.25 * (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:])
            if y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1 and _is_numeric(y):
                y = 0.25 * (y[:-1, :-1] + y[:-1, 1:] + y[1:, :-1] + y[1:, 1:])
    elif z.shape[-1] != xlen or z.shape[0] != ylen:
        # Helpful error message
        raise ValueError(
            f"Input shapes x {x.shape} and y {y.shape} "
            f"must match z centers {z.shape} "
            f"or z borders {tuple(i+1 for i in z.shape)}."
        )
    return x, y


# Input argument processing
def _from_data(data, *args):
    """
    Try to convert positional `key` arguments to `data[key]`. If argument is string
    it could be a valid positional argument like `fmt` so do not raise error.
    """
    if data is None:
        return
    args = list(args)
    for i, arg in enumerate(args):
        if isinstance(arg, str):
            try:
                array = data[arg]
            except KeyError:
                pass
            else:
                args[i] = array
    return args


def _preprocess_or_redirect(*keys, keywords=None, allow_extra=True):
    """
    Redirect internal plotting calls to native matplotlib methods. Also convert
    keyword args to positional and pass arguments through 'data' dictionary.
    """
    # Keyword arguments processed through 'data'
    # Positional arguments are always processed through data
    keywords = keywords or ()
    if isinstance(keywords, str):
        keywords = (keywords,)

    def _decorator(func):
        name = func.__name__
        from . import _kwargs_to_args

        @functools.wraps(func)
        def _preprocess_or_redirect(self, *args, **kwargs):
            if getattr(self, "_internal_call", None):
                # Redirect internal matplotlib call to native function
                from ..axes import PlotAxes

                func_native = getattr(super(PlotAxes, self), name)
                return func_native(*args, **kwargs)
            else:
                # Impose default coordinate system
                from ..constructor import Proj

                if self._name == "basemap" and name in BASEMAP_FUNCS:
                    if kwargs.get("latlon", None) is None:
                        kwargs["latlon"] = True
                if self._name == "cartopy" and name in CARTOPY_FUNCS:
                    if kwargs.get("transform", None) is None:
                        kwargs["transform"] = PlateCarree()
                    else:
                        kwargs["transform"] = Proj(kwargs["transform"])

                # Process data args
                # NOTE: Raises error if there are more args than keys
                args, kwargs = _kwargs_to_args(
                    keys, *args, allow_extra=allow_extra, **kwargs
                )
                data = kwargs.pop("data", None)
                if data is not None:
                    args = _from_data(data, *args)
                    for key in set(keywords) & set(kwargs):
                        kwargs[key] = _from_data(data, kwargs[key])

                # Auto-setup matplotlib with the input unit registry
                _load_objects()
                for arg in args:
                    if ndarray is not DataArray and isinstance(arg, DataArray):
                        arg = arg.data
                    if ndarray is not Quantity and isinstance(arg, Quantity):
                        ureg = getattr(arg, "_REGISTRY", None)
                        if hasattr(ureg, "setup_matplotlib"):
                            ureg.setup_matplotlib(True)

                # Sanitize colors
                # Lookup can be a 2-tuple with colormap, and integer denoting the color from that map
                # it also sanitizes if the second color is a float. Note that for float, the color
                # will be picked proportional to idx/255
                from .. import colors as pcolor

                color = kwargs.pop("color", None)
                if isinstance(color, tuple) and len(color) == 2:
                    cmap, color = color
                    color = pcolor._cmap_database.get_cmap(cmap)(color)
                if color is not None:
                    kwargs["color"] = color
                # Call main function
                return func(self, *args, **kwargs)  # call unbound method

        return _preprocess_or_redirect

    return _decorator


# Stats utiltiies
def _dist_clean(distribution):
    """
    Clean the distribution data for processing by `boxplot` or `violinplot`.
    Handles np.ndarrays where the ndarray is a list of lists of variable sizes.
    """
    if isinstance(distribution, np.ndarray):
        if distribution.dtype == object:
            # Handle list of lists with variable sizes
            return tuple(
                np.array(sublist, dtype=float)
                for sublist in distribution
                if len(sublist) > 0
            )
        else:
            # Handle regular numpy arrays
            if distribution.ndim == 1:
                distribution = distribution[:, None]
            distribution, units = _to_masked_array(distribution)  # no copy needed
            distribution = tuple(
                distribution[..., i].compressed() for i in range(distribution.shape[-1])
            )
            if units is not None:
                distribution = tuple(dist * units for dist in distribution)
            return distribution
    elif isinstance(distribution, list):
        # Handle list of lists directly
        return tuple(
            np.array(sublist, dtype=float)
            for sublist in distribution
            if len(sublist) > 0
        )
    else:
        raise ValueError("Input must be a numpy array or a list of lists")


def _dist_reduce(data, *, mean=None, means=None, median=None, medians=None, **kwargs):
    """
    Reduce statistical distributions to means and medians. Tack on a
    distribution keyword argument for processing down the line.
    """
    # TODO: Permit 3D array with error dimension coming first
    means = _not_none(mean=mean, means=means)
    medians = _not_none(median=median, medians=medians)
    if means and medians:
        warnings._warn_proplot(
            "Cannot have both means=True and medians=True. Using former."
        )
        medians = None
    if means or medians:
        distribution, units = _to_masked_array(data)
        distribution = distribution.filled()
        if distribution.ndim != 2:
            raise ValueError(
                f"Expected 2D array for means=True. Got {distribution.ndim}D."
            )
        if units is not None:
            distribution = distribution * units
        if means:
            data = np.nanmean(distribution, axis=0)
        else:
            data = np.nanmedian(distribution, axis=0)
        kwargs["distribution"] = distribution

    # Save argument passed to _error_bars
    return (data, kwargs)


def _dist_range(
    data,
    distribution,
    *,
    errdata=None,
    absolute=False,
    label=False,
    stds=None,
    pctiles=None,
    stds_default=None,
    pctiles_default=None,
):
    """
    Return a plottable characteristic range for the statistical distribution
    relative to the input coordinate (generally a mean or median).
    """
    # Parse stds arguments
    # NOTE: Have to guard against "truth value of an array is ambiguous" errors
    if stds is True:
        stds = stds_default
    elif stds is False or stds is None:
        stds = None
    else:
        stds = np.atleast_1d(stds)
        if stds.size == 1:
            stds = sorted((-stds.item(), stds.item()))
        elif stds.size != 2:
            raise ValueError("Expected scalar or length-2 stdev specification.")

    # Parse pctiles arguments
    if pctiles is True:
        pctiles = pctiles_default
    elif pctiles is False or pctiles is None:
        pctiles = None
    else:
        pctiles = np.atleast_1d(pctiles)
        if pctiles.size == 1:
            delta = (100 - pctiles.item()) / 2.0
            pctiles = sorted((delta, 100 - delta))
        elif pctiles.size != 2:
            raise ValueError("Expected scalar or length-2 pctiles specification.")

    # Incompatible settings
    if distribution is None and any(_ is not None for _ in (stds, pctiles)):
        raise ValueError(
            "To automatically compute standard deviations or percentiles on "
            "columns of data you must pass means=True or medians=True."
        )
    if stds is not None and pctiles is not None:
        warnings._warn_proplot(
            "Got both a standard deviation range and a percentile range for "
            "auto error indicators. Using the standard deviation range."
        )
        pctiles = None
    if distribution is not None and errdata is not None:
        stds = pctiles = None
        warnings._warn_proplot(
            "You explicitly provided the error bounds but also requested "
            "automatically calculating means or medians on data columns. "
            'It may make more sense to use the "stds" or "pctiles" keyword args '
            "and have *proplot* calculate the error bounds."
        )

    # Compute error data in format that can be passed to maxes.Axes.errorbar()
    # NOTE: Include option to pass symmetric deviation from central points
    if errdata is not None:
        # Manual error data
        if data.ndim != 1:
            raise ValueError(
                "Passing both 2D data coordinates and 'errdata' is not yet supported."
            )
        label_default = "uncertainty"
        err = _to_numpy_array(errdata)
        if (
            err.ndim not in (1, 2)
            or err.shape[-1] != data.size
            or err.ndim == 2
            and err.shape[0] != 2
        ):
            raise ValueError(
                f"Input 'errdata' has shape {err.shape}. Expected (2, {data.size})."
            )
        if err.ndim == 1:
            abserr = err
            err = np.empty((2, err.size))
            err[0, :] = data - abserr  # translated back to absolute deviations below
            err[1, :] = data + abserr
    elif stds is not None:
        # Standard deviations
        # NOTE: Invalid values were handled by _dist_reduce
        label_default = rf"{abs(stds[1])}$\sigma$ range"
        stds = _to_numpy_array(stds)[:, None]
        err = data + stds * np.nanstd(distribution, axis=0)
    elif pctiles is not None:
        # Percentiles
        # NOTE: Invalid values were handled by _dist_reduce
        label_default = f"{pctiles[1] - pctiles[0]}% range"
        err = np.nanpercentile(distribution, pctiles, axis=0)
    else:
        warnings._warn_proplot(
            "Error indications are missing from the dataset reduced by a "
            "mean or median operation. Consider passing e.g. bars=True."
        )
        err = None

    # Adjust error bounds
    if err is not None and not absolute:  # for errorbar() ingestion
        err = err - data
        err[0, :] *= -1  # absolute deviations from central points

    # Apply legend entry
    if isinstance(label, str):
        pass
    elif label:  # e.g. label=True says to use a default label
        label = label_default
    else:
        label = None

    return err, label


def _safe_mask(mask, *args):
    """
    Safely apply the mask to the input arrays, accounting for existing masked
    or invalid values. Values matching ``False`` are set to `np.nan`.
    """
    # NOTE: Could also convert unmasked data to masked. But other way around is
    # easier becase np.ma gives us correct fill values for data subtypes.
    _load_objects()
    invalid = ~mask  # True if invalid
    args_masked = []
    for data in args:
        data, units = _to_masked_array(data, copy=True)
        nan = data.fill_value
        data = data.filled()
        if data.size > 1 and data.shape != invalid.shape:
            raise ValueError(
                f"Mask shape {mask.shape} incompatible with array shape {data.shape}."
            )
        if data.size == 1 or invalid.size == 1:  # NOTE: happens with _restrict_inbounds
            pass
        elif invalid.size == 1:
            data = nan if invalid.item() else data
        elif data.size > 1:
            data[invalid] = nan
        if units is not None:
            data = data * units
        args_masked.append(data)
    return args_masked[0] if len(args_masked) == 1 else args_masked


def _safe_range(data, lo=0, hi=100):
    """
    Safely return the minimum and maximum (default) or percentile range accounting
    for masked values. Use min and max functions when possible for speed. Return
    ``None`` if we fail to get a valid range.
    """
    _load_objects()
    data, units = _to_masked_array(data)
    data = data.compressed()  # remove all invalid values
    min_ = max_ = None
    if data.size:
        min_ = np.min(data) if lo <= 0 else np.percentile(data, lo)
        if hasattr(min_, "dtype") and np.issubdtype(min_.dtype, np.integer):
            min_ = np.float64(min_)
        try:
            is_finite = np.isfinite(min_)
        except TypeError:
            is_finite = True
        if not is_finite:
            min_ = None
        elif units is not None:
            min_ *= units
    if data.size:
        max_ = np.max(data) if hi >= 100 else np.percentile(data, hi)
        if hasattr(max_, "dtype") and np.issubdtype(max_.dtype, np.integer):
            max_ = np.float64(max_)
        try:
            is_finite = np.isfinite(min_)
        except TypeError:
            is_finite = True
        if not is_finite:
            max_ = None
        elif units is not None:
            max_ *= units
    return min_, max_


# Metadata utilities
def _meta_coords(*args, which="x", **kwargs):
    """
    Return the index arrays associated with string coordinates and
    keyword arguments updated with index locators and formatters.
    """
    # NOTE: Why FixedLocator and not IndexLocator? The ticks chosen by the latter
    # depend on other plotted content.
    # NOTE: Why IndexFormatter and not FixedFormatter? The former ensures labels
    # correspond to indices while the latter can mysteriously truncate labels.
    from ..constructor import Formatter, Locator

    res = []
    for data in args:
        data = _to_duck_array(data)
        if not _is_categorical(data):
            res.append(data)
            continue
        if data.ndim > 1:
            raise ValueError("Non-1D string coordinate input is unsupported.")
        ticks = np.arange(len(data))
        labels = list(map(str, data))
        kwargs.setdefault(which + "locator", Locator(ticks))
        kwargs.setdefault(which + "formatter", Formatter(labels, index=True))
        kwargs.setdefault(which + "minorlocator", Locator("null"))
        res.append(ticks)  # use these as data coordinates
    return (*res, kwargs)


def _meta_labels(data, axis=0, always=True):
    """
    Return the array-like "labels" along axis `axis`. If `always` is ``False``
    we return ``None`` for simple ndarray input.
    """
    # NOTE: Previously inferred 'axis 1' metadata of 1D variable using the
    # data values metadata but that is incorrect. The paradigm for 1D plots
    # is we have row coordinates representing x, data values representing y,
    # and column coordinates representing individual series.
    _load_objects()
    labels = None
    if axis not in (0, 1, 2):
        raise ValueError(f"Invalid axis {axis}.")
    if isinstance(data, (ndarray, Quantity)):
        if not always:
            pass
        elif axis < data.ndim:
            labels = np.arange(data.shape[axis])
        else:  # requesting 'axis 1' on a 1D array
            labels = np.array([0])
    # Xarray object
    # NOTE: Even if coords not present .coords[dim] auto-generates indices
    elif isinstance(data, DataArray):
        if axis < data.ndim:
            labels = data.coords[data.dims[axis]]
        elif not always:
            pass
        else:
            labels = np.array([0])
    # Pandas object
    elif isinstance(data, (DataFrame, Series, Index)):
        if axis == 0 and isinstance(data, (DataFrame, Series)):
            labels = data.index
        elif axis == 1 and isinstance(data, (DataFrame,)):
            labels = data.columns
        elif not always:
            pass
        else:  # beyond dimensionality
            labels = np.array([0])
    # Everything else
    # NOTE: Ensure data is at least 1D in _to_duck_array so this covers everything
    else:
        raise ValueError(f"Unrecognized array type {type(data)}.")
    return labels


def _meta_title(data, include_units=True):
    """
    Return the "title" of an array-like object with metadata.
    Include units in the title if `include_units` is ``True``.
    """
    _load_objects()
    title = units = None
    if isinstance(data, ndarray):
        pass
    # Xarray object with possible long_name, standard_name, and units attributes.
    # Output depends on if units is True. Prefer long_name (come last in loop).
    elif isinstance(data, DataArray):
        title = getattr(data, "name", None)
        for key in ("standard_name", "long_name"):
            title = data.attrs.get(key, title)
        if include_units:
            units = _meta_units(data)
    # Pandas object. DataFrame has no native name attribute but user can add one
    # See: https://github.com/pandas-dev/pandas/issues/447
    elif isinstance(data, (DataFrame, Series, Index)):
        title = getattr(data, "name", None) or None
    # Pint Quantity
    elif isinstance(data, Quantity):
        if include_units:
            units = _meta_units(data)
    # Add units or return units alone
    if title and units:
        title = f"{title} ({units})"
    else:
        title = title or units
    if title is not None:
        title = str(title).strip()
    return title


def _meta_units(data):
    """
    Get the unit string from the `xarray.DataArray` attributes or the
    `pint.Quantity`. Format the latter with :rcraw:`unitformat`.
    """
    _load_objects()
    # Get units from the attributes
    if ndarray is not DataArray and isinstance(data, DataArray):
        units = data.attrs.get("units", None)
        data = data.data
        if units is not None:
            return units
    # Get units from the quantity
    if ndarray is not Quantity and isinstance(data, Quantity):
        from ..config import rc

        fmt = rc.unitformat
        try:
            units = format(data.units, fmt)
        except (TypeError, ValueError):
            warnings._warn_proplot(
                f"Failed to format pint quantity with format string {fmt!r}."
            )
        else:
            if "L" in fmt:  # auto-apply LaTeX math indicator
                units = "$" + units + "$"
            return units


# Geographic utiltiies
def _geo_basemap_1d(x, *ys, xmin=-180, xmax=180):
    """
    Fix basemap geographic 1D data arrays.
    """
    ys = _geo_clip(*ys)
    x_orig, ys_orig, ys = x, ys, []
    for y_orig in ys_orig:
        x, y = _geo_inbounds(x_orig, y_orig, xmin=xmin, xmax=xmax)
        ys.append(y)
    return (x, *ys)


def _geo_basemap_2d(x, y, *zs, xmin=-180, xmax=180, globe=False):
    """
    Fix basemap geographic 2D data arrays.
    """
    y = _geo_clip(y)
    x_orig, y_orig, zs_orig, zs = x, y, zs, []
    for z_orig in zs_orig:
        x, y, z = x_orig, y_orig, z_orig
        x, z = _geo_inbounds(x, z, xmin=xmin, xmax=xmax)
        if globe and z is not None and x.ndim == 1 and y.ndim == 1:
            x, y, z = _geo_globe(x, y, z, xmin=xmin, modulo=False)
        zs.append(z)
    return (x, y, *zs)


def _geo_cartopy_1d(x, *ys):
    """
    Fix cartopy geographic 1D data arrays.
    """
    ys = _geo_clip(ys)
    return (x, *ys)


def _geo_cartopy_2d(x, y, *zs, globe=False):
    """
    Fix cartopy geographic 2D data arrays.
    """
    y = _geo_clip(y)
    x_orig, y_orig, zs_orig = x, y, zs
    zs = []
    for z_orig in zs_orig:
        x, y, z = x_orig, y_orig, z_orig
        if globe and z is not None and x.ndim == 1 and y.ndim == 1:
            x, y, z = _geo_globe(x, y, z, modulo=True)
        zs.append(z)
    return (x, y, *zs)


def _geo_clip(*ys):
    """
    Ensure latitudes fall within ``-90`` to ``90``. Important if we
    add graticule edges with `edges`.
    """
    ys = tuple(np.clip(y, -90, 90) for y in ys)
    return ys[0] if len(ys) == 1 else ys


def _geo_inbounds(x, y, xmin=-180, xmax=180):
    """
    Fix conflicts with map coordinates by rolling the data to fall between the
    minimum and maximum longitudes and masking out-of-bounds data points.
    """
    # Roll in same direction if some points on right-edge extend
    # more than 360 above min longitude; *they* should be on left side
    if x.ndim != 1:
        return x, y
    lonroll = np.where(x > xmin + 360)[0]  # tuple of ids
    if lonroll.size:  # non-empty
        roll = x.size - lonroll.min()
        x = np.roll(x, roll)
        y = np.roll(y, roll, axis=-1)
        x[:roll] -= 360  # make monotonic
    # Set NaN where data not in range xmin, xmax. Must be done for regional smaller
    # projections or get weird side-effects from valid data outside boundaries
    y, units = _to_masked_array(y)
    nan = y.fill_value
    y = y.filled()
    if not y.shape:
        pass
    elif x.size - 1 == y.shape[-1]:  # test western/eastern grid cell edges
        mask = (x[1:] < xmin) | (x[:-1] > xmax)
        y[..., mask] = nan
    elif x.size == y.shape[-1]:  # test the centers and pad by one for safety
        (where,) = np.where((x < xmin) | (x > xmax))
        y[..., where[1:-1]] = nan
    return x, y


def _geo_globe(x, y, z, xmin=-180, modulo=False):
    """
    Ensure global coverage by fixing gaps over poles and across
    longitude seams. Increases the size of the arrays.
    """
    # Cover gaps over poles by appending polar data
    with np.errstate(all="ignore"):
        p1 = np.mean(z[0, :])  # do not ignore NaN if present
        p2 = np.mean(z[-1, :])
    ps = (-90, 90) if (y[0] < y[-1]) else (90, -90)
    z1 = np.repeat(p1, z.shape[1])
    z2 = np.repeat(p2, z.shape[1])
    y = ma.concatenate((ps[:1], y, ps[1:]))
    z = ma.concatenate((z1[None, :], z, z2[None, :]), axis=0)
    # Cover gaps over cartopy longitude seam
    # Ensure coordinates span 360 after modulus
    if modulo:
        if x[0] % 360 != (x[-1] + 360) % 360:
            x = ma.concatenate((x, (x[0] + 360,)))
            z = ma.concatenate((z, z[:, :1]), axis=1)
    # Cover gaps over basemap longitude seam
    # Ensure coordinates span exactly 360
    else:
        # Interpolate coordinate centers to seam. Size possibly augmented by 2
        if x.size == z.shape[1]:
            if x[0] + 360 != x[-1]:
                xi = np.array([x[-1], x[0] + 360])  # input coordinates
                xq = xmin + 360  # query coordinate
                zq = ma.concatenate((z[:, -1:], z[:, :1]), axis=1)
                zq = (zq[:, :1] * (xi[1] - xq) + zq[:, 1:] * (xq - xi[0])) / (
                    xi[1] - xi[0]
                )  # noqa: E501
                x = ma.concatenate(((xmin,), x, (xmin + 360,)))
                z = ma.concatenate((zq, z, zq), axis=1)
        # Extend coordinate edges to seam. Size possibly augmented by 1.
        elif x.size - 1 == z.shape[1]:
            if x[0] != xmin:
                x = ma.append(xmin, x)
                x[-1] = xmin + 360
                z = ma.concatenate((z[:, -1:], z), axis=1)
        else:
            raise ValueError("Unexpected shapes of coordinates or data arrays.")
    return x, y, z
