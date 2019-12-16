#!/usr/bin/env python3
# Import everything into the top-level module namespace
# Monkey patch warnings format for warnings issued by ProPlot, make sure to
# detect if this is just a matplotlib warning traced back to ProPlot code by
# testing whether the warned line contains "warnings.warn"
# See: https://stackoverflow.com/a/2187390/4970632
# For internal warning call signature:
# https://docs.python.org/3/library/warnings.html#warnings.showwarning
# For default warning source code see:
# https://github.com/python/cpython/blob/master/Lib/warnings.py
import os as _os
import pkg_resources as _pkg
from .utils import _benchmark
from .utils import *  # noqa: F401 F403
with _benchmark('total time'):
    with _benchmark('styletools'):
        from .styletools import *  # noqa: F401 F403
    with _benchmark('rctools'):
        from .rctools import *  # noqa: F401 F403
    with _benchmark('axistools'):
        from .axistools import *  # noqa: F401 F403
    with _benchmark('wrappers'):
        from .wrappers import *  # noqa: F401 F403
    with _benchmark('projs'):
        from .projs import *  # noqa: F401 F403
    with _benchmark('axes'):
        from .axes import *  # noqa: F401 F403
    with _benchmark('subplots'):
        from .subplots import *  # noqa: F401 F403


# Initialize customization folders
_rc_folder = _os.path.join(_os.path.expanduser('~'), '.proplot')
if not _os.path.isdir(_rc_folder):
    _os.mkdir(_rc_folder)
for _rc_sub in ('cmaps', 'cycles', 'colors', 'fonts'):
    _rc_sub = _os.path.join(_rc_folder, _rc_sub)
    if not _os.path.isdir(_rc_sub):
        _os.mkdir(_rc_sub)

# SCM versioning
name = 'proplot'
try:
    version = __version__ = _pkg.get_distribution(__name__).version
except _pkg.DistributionNotFound:
    version = __version__ = 'unknown'
