#!/usr/bin/env python3
"""
Utilties related to legends and colorbars.
"""
import matplotlib.artist as martist
import matplotlib.colorbar as mcolorbar
import matplotlib.legend as mlegend  # noqa: F401
import matplotlib.ticker as mticker
import numpy as np

from . import ic  # noqa: F401
from . import warnings

# Global constants
REMOVE_AFTER_FLUSH = (
    "pad",
    "space",
    "width",
    "length",
    "shrink",
    "align",
    "queue",
)
GUIDE_ALIASES = (
    ("title", "label"),
    ("locator", "ticks"),
    ("format", "formatter", "ticklabels"),
)


def _add_guide_kw(name, kwargs, **opts):
    """
    Add to the `colorbar_kw` or `legend_kw` dict if there are no conflicts.
    """
    # NOTE: Here we *do not* want to overwrite properties in dictionary. Indicates
    # e.g. default locator inferred from levels or default title inferred from metadata.
    attr = f"{name}_kw"
    if not opts:
        return
    if not kwargs.get(attr, None):
        kwargs[attr] = {}  # permit e.g. colorbar_kw=None
    guide_kw = kwargs[attr]
    _update_kw(guide_kw, overwrite=False, **opts)


def _cache_guide_kw(obj, name, kwargs):
    """
    Cache settings on the object from the input keyword arguments.
    """
    # NOTE: Here we overwrite the hidden dictionary if it already exists.
    # This is only called once in _update_guide() so its fine.
    try:
        setattr(obj, f"_{name}_kw", kwargs)
    except AttributeError:
        pass
    if isinstance(obj, (tuple, list, np.ndarray)):
        for member in obj:
            _cache_guide_kw(member, name, kwargs)


def _flush_guide_kw(obj, name, kwargs):
    """
    Flux settings cached on the object into the keyword arguments.
    """
    # NOTE: Here we *do not* overwrite properties in the dictionary by default.
    # Indicates e.g. calling colorbar(extend='both') after pcolor(extend='neither').
    # NOTE: Previously had problems reusing same keyword arguments for more than one
    # colorbar() because locator or formatter axis would get reset. Old solution was
    # to delete the _guide_kw but that destroyed default behavior. New solution is
    # to keep _guide_kw but have constructor functions return shallow copies.
    guide_kw = getattr(obj, f"_{name}_kw", None)
    if guide_kw:
        _update_kw(kwargs, overwrite=False, **guide_kw)
        for key in REMOVE_AFTER_FLUSH:
            guide_kw.pop(key, None)
    if isinstance(obj, (tuple, list, np.ndarray)):
        for member in obj:  # possibly iterate over matplotlib tuple/list subclasses
            _flush_guide_kw(member, name, kwargs)
    return kwargs


def _update_kw(kwargs, overwrite=False, **opts):
    """
    Add the keyword arguments to the dictionary if not already present.
    """
    for key, value in opts.items():
        if value is None:
            continue
        keys = tuple(k for opts in GUIDE_ALIASES for k in opts if key in opts)
        keys = keys or (key,)  # e.g. 'extend' or something
        keys_found = tuple(key for key in keys if kwargs.get(key) is not None)
        if not keys_found:
            kwargs[key] = value
        elif overwrite:  # overwrite existing key
            kwargs[keys_found[0]] = value


def _iter_children(*args):
    """
    Iterate through `_children` of `HPacker`, `VPacker`, and `DrawingArea`.
    This is used to update legend handle properties.
    """
    for arg in args:
        if hasattr(arg, "_children"):
            yield from _iter_children(*arg._children)
        elif arg is not None:
            yield arg


def _iter_iterables(*args):
    """
    Iterate over arbitrary nested lists of iterables. Used for deciphering legend input.
    Things can get complicated with e.g. bar colunns plus negative-positive colors.
    """
    for arg in args:
        if np.iterable(arg):
            yield from _iter_iterables(*arg)
        elif arg is not None:
            yield arg


def _update_ticks(self, manual_only=False):
    """
    Refined colorbar tick updater without subclassing.
    """
    # TODO: Add this to generalized colorbar subclass?
    # NOTE: Matplotlib 3.5+ does not define _use_auto_colorbar_locator since
    # ticks are always automatically adjusted by its colorbar subclass. This
    # override is thus backwards and forwards compatible.
    attr = "_use_auto_colorbar_locator"
    if not hasattr(self, attr) or getattr(self, attr)():
        if manual_only:
            pass
        else:
            mcolorbar.Colorbar.update_ticks(self)  # AutoMinorLocator auto updates
    else:
        mcolorbar.Colorbar.update_ticks(self)  # update necessary
        minorlocator = getattr(self, "minorlocator", None)
        if minorlocator is None:
            pass
        elif hasattr(self, "_ticker"):
            ticks, *_ = self._ticker(self.minorlocator, mticker.NullFormatter())
            axis = self.ax.yaxis if self.orientation == "vertical" else self.ax.xaxis
            axis.set_ticks(ticks, minor=True)
            axis.set_ticklabels([], minor=True)
        else:
            warnings._warn_proplot(
                f"Cannot use user-input colorbar minor locator {minorlocator!r} (older matplotlib version). Turning on minor ticks instead."
            )  # noqa: E501
            self.minorlocator = None
            self.minorticks_on()  # at least turn them on


class _InsetColorbar(martist.Artist):
    """
    Legend-like class for managing inset colorbars.
    """

    # TODO: Write this!


class _CenteredLegend(martist.Artist):
    """
    Legend-like class for managing centered-row legends.
    """

    # TODO: Write this!
