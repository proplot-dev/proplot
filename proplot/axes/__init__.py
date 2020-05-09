#!/usr/bin/env python3
"""
The axes classes used for all ProPlot figures.
"""
from .base import Axes  # noqa: F401
from ..internals import warnings
from .cartesian import CartesianAxes
from .polar import PolarAxes
from .geo import GeoAxes, BasemapAxes, CartopyAxes  # noqa: F401
XYAxes = warnings._rename_obj('XYAxes', CartesianAxes)
ProjAxes = warnings._rename_obj('ProjAxes', GeoAxes)

import matplotlib.projections as mproj
mproj.register_projection(CartesianAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)

# Prevent importing module names and set order of appearance for objects
__all__ = [
    'Axes', 'CartesianAxes', 'PolarAxes',
    'GeoAxes', 'CartopyAxes', 'BasemapAxes',
    'ProjAxes', 'XYAxes',  # deprecated
]
