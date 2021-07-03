#!/usr/bin/env python3
"""
The "3D" axes class.
"""
from ..config import rc
from . import base

try:
    import mpl_toolkits.mplot3d as three
except ImportError:
    three = Axes3D = None

if three is not None:
    class Axes3D(base.Axes, three.Axes3D):
        """
        Simple mix-in of `~proplot.axes.Axes` with `~mpl_toolkits.mplot3d.Axes3D`.
        """
        #: The registered projection name.
        name = 'proplot_3d'

        def format(self, **kwargs):
            # No-op for now
            rc_kw, rc_mode, kwargs = self._parse_format(**kwargs)
            with rc.context(rc_kw, mode=rc_mode):
                return super().format(**kwargs)
