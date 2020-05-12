#!/usr/bin/env python3
"""
The axes classes used for all ProPlot figures.
"""
from . import plot
from .plot import *  # noqa: F401, F403
from .base import Axes  # noqa: F401
from .cartesian import CartesianAxes
from .polar import PolarAxes
from .geo import GeoAxes  # noqa: F401
from .geo import BasemapAxes, CartopyAxes
from ..internals import warnings
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
__all__.extend(plot.__all__)
