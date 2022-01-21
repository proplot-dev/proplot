#!/usr/bin/env python3
"""
A succinct matplotlib wrapper for making beautiful, publication-quality graphics.
"""
# SCM versioning
import pkg_resources as pkg
name = 'proplot'
try:
    version = __version__ = pkg.get_distribution(__name__).version
except pkg.DistributionNotFound:
    version = __version__ = 'unknown'

# Import optional dependencies now to isolate import times
from .internals.benchmarks import _benchmark
with _benchmark('pyplot'):
    from matplotlib import pyplot  # noqa: F401
with _benchmark('cartopy'):
    try:
        import cartopy  # noqa: F401
    except ImportError:
        pass
with _benchmark('basemap'):
    try:
        from mpl_toolkits import basemap  # noqa: F401
    except ImportError:
        pass

# Import everything to top level
with _benchmark('config'):
    from .config import *  # noqa: F401 F403
with _benchmark('crs'):
    from .crs import *  # noqa: F401 F403
with _benchmark('utils'):
    from .utils import *  # noqa: F401 F403
with _benchmark('colors'):
    from .colors import *  # noqa: F401 F403
with _benchmark('ticker'):
    from .ticker import *  # noqa: F401 F403
with _benchmark('scale'):
    from .scale import *  # noqa: F401 F403
with _benchmark('axes'):
    from .axes import *  # noqa: F401 F403
with _benchmark('gridspec'):
    from .gridspec import *  # noqa: F401 F403
with _benchmark('figure'):
    from .figure import *  # noqa: F401 F403
with _benchmark('constructor'):
    from .constructor import *  # noqa: F401 F403
with _benchmark('ui'):
    from .ui import *  # noqa: F401 F403
with _benchmark('demos'):
    from .demos import *  # noqa: F401 F403
    from .artist import *  # noqa: F401 F403

# Dynamically add registered classes to top-level namespace
from .constructor import NORMS, LOCATORS, FORMATTERS, SCALES, PROJS
_globals = globals()
for _src in (NORMS, LOCATORS, FORMATTERS, SCALES, PROJS):
    for _key, _cls in _src.items():
        if isinstance(_cls, type):  # i.e. not a scale preset
            _globals[_cls.__name__] = _cls  # may overwrite proplot names

# Register objects
from .config import register_cmaps, register_cycles, register_colors, register_fonts
with _benchmark('cmaps'):
    register_cmaps(default=True)
with _benchmark('cycles'):
    register_cycles(default=True)
with _benchmark('colors'):
    register_colors(default=True)
with _benchmark('fonts'):
    register_fonts(default=True)

# Validate colormap names and propagate 'cycle' to 'axes.prop_cycle'
# NOTE: cmap.sequential also updates siblings 'cmap' and 'image.cmap'
from .config import rc
from .internals import rcsetup, warnings
rcsetup.VALIDATE_REGISTERED_CMAPS = True
for _key in ('cycle', 'cmap.sequential', 'cmap.diverging', 'cmap.cyclic', 'cmap.qualitative'):  # noqa: E501
    try:
        rc[_key] = rc[_key]
    except ValueError as err:
        warnings._warn_proplot(f'Invalid user rc file setting: {err}')
        rc[_key] = 'Greys'  # fill value

# Validate color names now that colors are registered
# NOTE: This updates all settings with 'color' in name (harmless if it's not a color)
from .config import rc_proplot, rc_matplotlib
rcsetup.VALIDATE_REGISTERED_COLORS = True
for _src in (rc_proplot, rc_matplotlib):
    for _key in _src:  # loop through unsynced properties
        if 'color' not in _key:
            continue
        try:
            _src[_key] = _src[_key]
        except ValueError as err:
            warnings._warn_proplot(f'Invalid user rc file setting: {err}')
            _src[_key] = 'black'  # fill value
