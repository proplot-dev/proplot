#!/usr/bin/env python3
"""
The axes classes used for all ProPlot figures.
"""
import matplotlib.projections as mproj

from ..internals import warnings
from . import plot
from .base import Axes  # noqa: F401
from .cartesian import CartesianAxes
from .geo import GeoAxes  # noqa: F401
from .geo import BasemapAxes, CartopyAxes
from .plot import *  # noqa: F401, F403
from .polar import PolarAxes
from .three import Axes3D  # noqa: F401

XYAxes, ProjAxes = warnings._rename_objs(
    '0.6',
    XYAxes=CartesianAxes,
    ProjAxes=GeoAxes,
)

mproj.register_projection(CartesianAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)
if Axes3D is not None:
    mproj.register_projection(Axes3D)

# Prevent importing module names and set order of appearance for objects
__all__ = [
    'Axes', 'CartesianAxes', 'PolarAxes',
    'GeoAxes', 'CartopyAxes', 'BasemapAxes',
    'ProjAxes', 'XYAxes',  # deprecated
]
__all__.extend(plot.__all__)  # document wrappers as part of proplot/axes submodule
