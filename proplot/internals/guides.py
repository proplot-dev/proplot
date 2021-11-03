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


def _fill_guide_kw(kwargs, **pairs):
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
        keys = tuple(a for group in aliases for a in group if key in group)  # may be ()
        if not any(kwargs.get(key) is not None for key in keys):  # note any(()) is True
            kwargs[key] = value


def _guide_kw_from_obj(obj, name, kwargs):
    """
    Add to the dict from settings stored on the object if there are no conflicts.
    """
    pairs = getattr(obj, f'_{name}_kw', None)
    pairs = pairs or {}  # needed for some reason
    _fill_guide_kw(kwargs, **pairs)
    if isinstance(obj, (tuple, list, np.ndarray)):
        for iobj in obj:  # possibly iterate over matplotlib tuple/list subclasses
            _guide_kw_from_obj(iobj, name, kwargs)
    return kwargs


def _guide_kw_to_obj(obj, name, kwargs):
    """
    Add the guide keyword dict to the objects.
    """
    try:
        setattr(obj, f'_{name}_kw', kwargs)
    except AttributeError:
        pass
    if isinstance(obj, (tuple, list, np.ndarray)):
        for iobj in obj:
            _guide_kw_to_obj(iobj, name, kwargs)


def _guide_kw_to_arg(name, kwargs, **pairs):
    """
    Add to the `colorbar_kw` or `legend_kw` dict if there are no conflicts.
    """
    kw = kwargs.setdefault(f'{name}_kw', {})
    _fill_guide_kw(kw, **pairs)


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
    # WARNING: Important to guard against colorbar private API changes here
    use_auto = getattr(self, '_use_auto_colorbar_locator', lambda: True)
    if not use_auto():
        mcolorbar.Colorbar.update_ticks(self)  # update necessary
        minorlocator = getattr(self, 'minorlocator', None)
        if minorlocator is not None and hasattr(self, '_ticker'):
            ticks, *_ = self._ticker(self.minorlocator, mticker.NullFormatter())
            axis = self.ax.yaxis if self.orientation == 'vertical' else self.ax.xaxis
            axis.set_ticks(ticks, minor=True)
            axis.set_ticklabels([], minor=True)
    elif not manual_only:
        mcolorbar.Colorbar.update_ticks(self)  # here AutoMinorLocator auto updates


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
