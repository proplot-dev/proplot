#!/usr/bin/env python3
"""
The "3D" axes class.
"""
from . import base, shared

try:
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    Axes3D = object


class ThreeAxes(shared._SharedAxes, base.Axes, Axes3D):
    """
    Simple mix-in of `proplot.axes.Axes` with `~mpl_toolkits.mplot3d.axes3d.Axes3D`.

    Important
    ---------
    This axes subclass does not implement the `~proplot.axes.PlotAxes` plotting
    overrides. It can be used by passing ``proj='3d'`` or ``proj='three'``
    to axes-creation commands like `~proplot.figure.Figure.add_axes`,
    `~proplot.figure.Figure.add_subplot`, and `~proplot.figure.Figure.subplots`.
    """
    _name = 'three'
    _name_aliases = ('3d',)

    def __init__(self, *args, **kwargs):
        import mpl_toolkits.mplot3d  # noqa: F401 verify package is available
        kwargs.setdefault('alpha', 0.0)
        super().__init__(*args, **kwargs)
