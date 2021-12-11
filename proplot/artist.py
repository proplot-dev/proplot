#!/usr/bin/env python3
"""
Add dot-notation properties for matplotlib setter and getter functions.
"""
import warnings

import matplotlib.artist as martist
from matplotlib import MatplotlibDeprecationWarning

__all__ = []


PROPS_IGNORE = (
    # Axes props
    'axes',
    'figure',
    'xaxis',  # axes prop
    'yaxis',  # axes prop
    'zaxis',  # axes prop
    'units',  # axis prop
    'gridlines',
    # Method conflicts
    'legend',
    'tight_layout',
    # Class-level conflicts
    'contains',
    'zorder',
    'pickradius',  # line property related to 'contains'
    # Instance-level artist conflicts
    'label',
    'label_position',
    # Instance-level axes conflicts
    'cmap',
    'norm',
    'lines',
    'images',
    'title',  # TODO: use internal title handling
)


def _iter_subclasses(cls):
    """
    Iterate through all subclasses.
    """
    yield cls
    try:
        for subclass in cls.__subclasses__():
            yield from _iter_subclasses(subclass)
    except TypeError:
        pass


def _add_properties(cls):
    """
    Generate property definitions for every artist getter.
    """
    for attr in dir(cls):
        try:
            getter = getattr(cls, attr)
        except MatplotlibDeprecationWarning:
            continue
        if not callable(getter) or attr[:4] != 'get_':
            continue
        prop = attr[4:]
        if prop in PROPS_IGNORE:
            continue
        if hasattr(cls, prop):
            value = getattr(cls, prop)
            if not isinstance(value, property):  # i.e. this is not child of a class
                warnings._warn_proplot(f'Skipping property {prop!r}. Already exists as attribute.')  # noqa: E501
            continue
        args = [getter]  # property() function args
        setter = getattr(cls, 'set_' + prop, None)
        if callable(setter):
            args.append(setter)
        obj = property(*args, doc=getter.__doc__)
        setattr(cls, prop, obj)


# Apply properties
# NOTE: While we can guard against class-level attribute conflicts we *cannot* guard
# against instance-level attribute conflicts. Therefore this may never work.
for cls in _iter_subclasses(martist.Artist):
    with warnings.catch_warnings():
        warnings.simplefilter('error', MatplotlibDeprecationWarning)
        _add_properties(cls)
