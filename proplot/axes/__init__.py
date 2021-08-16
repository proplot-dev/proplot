#!/usr/bin/env python3
"""
The axes classes used for all ProPlot figures.
"""
import matplotlib.projections as mproj

from .base import Axes  # noqa: F401
from .cartesian import CartesianAxes
from .geo import GeoAxes  # noqa: F401
from .geo import _BasemapAxes, _CartopyAxes
from .plot import PlotAxes  # noqa: F401
from .polar import PolarAxes
from .shared import _SharedAxes  # noqa: F401
from .three import ThreeAxes  # noqa: F401

# Register projections
mproj.register_projection(CartesianAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(_BasemapAxes)
mproj.register_projection(_CartopyAxes)
mproj.register_projection(ThreeAxes)

# Prevent importing module names and set order of appearance for objects
__all__ = [
    'Axes',
    'PlotAxes',
    'CartesianAxes',
    'PolarAxes',
    'GeoAxes',
    'ThreeAxes',
]
