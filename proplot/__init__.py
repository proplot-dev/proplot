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
import warnings as _warnings
import pkg_resources as _pkg
from .utils import _benchmark
with _benchmark('total time'):
    from .utils import *  # noqa
    with _benchmark('styletools'):
        from .styletools import *  # noqa
    with _benchmark('rctools'):
        from .rctools import *  # noqa
    with _benchmark('axistools'):
        from .axistools import *  # noqa
    with _benchmark('wrappers'):
        from .wrappers import *  # noqa
    with _benchmark('projs'):
        from .projs import *  # noqa
    with _benchmark('axes'):
        from .axes import *  # noqa
    with _benchmark('subplots'):
        from .subplots import *  # noqa


def _warning_proplot(message, category, filename, lineno, line=None):
    """
    Format for warnings issued by ProPlot. If this is
    just a matplotlib warning traced back to ProPlot code the *default*
    warning format is used.
    See the `internal warning call signature \
<https://docs.python.org/3/library/warnings.html#warnings.showwarning>`__
    and the `default warning source code \
<https://github.com/python/cpython/blob/master/Lib/warnings.py>`__.
    """
    if line is None:
        try:
            import linecache
            line = linecache.getline(filename, lineno)
        except ModuleNotFoundError:
            pass
    if 'proplot' in filename and line is not None and 'warnings' in line:
        string = f'{filename}:{lineno}: ProPlotWarning: {message}'
    else:
        string = f'{filename}:{lineno}: {category.__name__}: {message}'
        if line is not None:
            string += ('\n' + line)  # default behavior
    return (string + '\n')  # must end in newline or not shown in IPython


# Apply monkeypatch
# See: https://stackoverflow.com/a/2187390/4970632
if _warnings.formatwarning is not _warning_proplot:
    _warnings.formatwarning = _warning_proplot


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
