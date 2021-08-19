#!/usr/bin/env python3
"""
The "3D" axes class.
"""
from . import plot, shared

try:
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    Axes3D = object


class ThreeAxes(shared._SharedAxes, plot.PlotAxes, Axes3D):
    """
    Simple mix-in of `proplot.axes.PlotAxes` with `~mpl_toolkits.mplot3d.axes3d.Axes3D`.
    """
    #: The registered projection name.
    name = 'proplot_three'

    def __init__(self, *args, **kwargs):
        import mpl_toolkits.mplot3d  # noqa: F401 verify package is available
        # Initialize axes
        super().__init__(*args, **kwargs)

    def _update_background(self, **kwargs):
        # Force the figure face color to the axes patch color or else the axes
        # look haphazardly thrown onto a square background patch and the spines
        # and labels bleed into the figure edge region.
        super()._update_background(**kwargs)
        self.figure.patch.set_facecolor(self.patch.get_facecolor())
