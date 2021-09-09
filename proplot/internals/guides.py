#!/usr/bin/env python3
"""
Utilties related to legends and colorbars.
"""
import numpy as np


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
