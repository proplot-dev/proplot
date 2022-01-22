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


def _fill_guide_kw(kwargs, overwrite=False, **pairs):
    """
    Add the keyword arguments to the dictionary if not already present.
    """
    aliases = (
        ('title', 'label'),
        ('locator', 'ticks'),
        ('format', 'formatter', 'ticklabels')
    )
    for key, value in pairs.items():
        if value is None:
            continue
        keys = tuple(k for opts in aliases for k in opts if key in opts)
        keys = keys or (key,)  # e.g. 'extend' or something
        keys_found = tuple(key for key in keys if kwargs.get(key) is not None)
        if not keys_found:
            kwargs[key] = value
        elif overwrite:  # overwrite existing key
            kwargs[keys_found[0]] = value


def _guide_kw_to_arg(name, kwargs, **pairs):
    """
    Add to the `colorbar_kw` or `legend_kw` dict if there are no conflicts.
    """
    # WARNING: Here we *do not* want to overwrite properties in the dictionary.
    # Indicates e.g. calling colorbar(extend='both') after pcolor(extend='neither').
    kw = kwargs.setdefault(f'{name}_kw', {})
    _fill_guide_kw(kw, overwrite=True, **pairs)


def _guide_kw_to_obj(obj, name, kwargs):
    """
    Store settings on the object from the input dict.
    """
    try:
        setattr(obj, f'_{name}_kw', kwargs)
    except AttributeError:
        pass
    if isinstance(obj, (tuple, list, np.ndarray)):
        for member in obj:
            _guide_kw_to_obj(member, name, kwargs)


def _guide_obj_to_kw(obj, name, kwargs):
    """
    Add to the dict from settings stored on the object if there are no conflicts.
    """
    # WARNING: Here we *do* want to overwrite properties in dictionary. Indicates
    # updating kwargs during parsing (probably only relevant for ax.parametric).
    pairs = getattr(obj, f'_{name}_kw', None)
    pairs = pairs or {}  # needed for some reason
    _fill_guide_kw(kwargs, overwrite=False, **pairs)
    if isinstance(obj, (tuple, list, np.ndarray)):
        for member in obj:  # possibly iterate over matplotlib tuple/list subclasses
            _guide_obj_to_kw(member, name, kwargs)
    return kwargs


def _iter_children(*args):
    """
    Iterate through `_children` of `HPacker`, `VPacker`, and `DrawingArea`.
    This is used to update legend handle properties.
    """
    for arg in args:
        if hasattr(arg, '_children'):
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
    use_auto = getattr(self, '_use_auto_colorbar_locator', lambda: True)
    if use_auto():
        if manual_only:
            pass
        else:
            mcolorbar.Colorbar.update_ticks(self)  # here AutoMinorLocator auto updates
    else:
        mcolorbar.Colorbar.update_ticks(self)  # update necessary
        minorlocator = getattr(self, 'minorlocator', None)
        if minorlocator is None:
            pass
        elif hasattr(self, '_ticker'):
            ticks, *_ = self._ticker(self.minorlocator, mticker.NullFormatter())
            axis = self.ax.yaxis if self.orientation == 'vertical' else self.ax.xaxis
            axis.set_ticks(ticks, minor=True)
            axis.set_ticklabels([], minor=True)
        else:
            warnings._warn_proplot(f'Cannot use user-input colorbar minor locator {minorlocator!r} (older matplotlib version). Turning on minor ticks instead.')  # noqa: E501
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
