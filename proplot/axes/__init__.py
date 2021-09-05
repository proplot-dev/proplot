#!/usr/bin/env python3
"""
The various axes classes used throughout ProPlot.
"""
import matplotlib.projections as mproj

from ..internals import context
from .base import Axes  # noqa: F401
from .cartesian import CartesianAxes
from .geo import GeoAxes  # noqa: F401
from .geo import _BasemapAxes, _CartopyAxes
from .plot import PlotAxes  # noqa: F401
from .polar import PolarAxes
from .shared import _SharedAxes  # noqa: F401
from .three import ThreeAxes  # noqa: F401

# Prevent importing module names and set order of appearance for objects
__all__ = [
    'Axes',
    'PlotAxes',
    'CartesianAxes',
    'PolarAxes',
    'GeoAxes',
    'ThreeAxes',
]

# Register projections with package prefix to avoid conflicts
# NOTE: We integrate with cartopy and basemap rather than using matplotlib's
# native projection system. Therefore axes names are not part of public API.
for _cls in (CartesianAxes, PolarAxes, _BasemapAxes, _CartopyAxes, ThreeAxes):
    for _name in (_cls._name, *_cls._name_aliases):
        with context._state_context(_cls, name='proplot_' + _name):
            mproj.register_projection(_cls)
