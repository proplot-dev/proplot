#!/usr/bin/env python3
"""
A powerful matplotlib wrapper for making beautiful, publication-quality graphics.
"""
# SCM versioning
import pkg_resources as pkg
name = 'proplot'
try:
    version = __version__ = pkg.get_distribution(__name__).version
except pkg.DistributionNotFound:
    version = __version__ = 'unknown'

# Import everything to top level
from .internals import benchmarks, docstring, rcsetup, warnings  # noqa: F401
with benchmarks._benchmark('imports'):
    from .config import *  # noqa: F401 F403
    from .crs import *  # noqa: F401 F403
    from .utils import *  # noqa: F401 F403
    from .colors import *  # noqa: F401 F403
    from .ticker import *  # noqa: F401 F403
    from .scale import *  # noqa: F401 F403
    from .gridspec import *  # noqa: F401 F403
    from .constructor import *  # noqa: F401 F403
    from .axes import *  # noqa: F401 F403
    from .figure import *  # noqa: F401 F403
    from .ui import *  # noqa: F401 F403
    from .demos import *  # noqa: F401 F403

# Dynamically add registered classes to top-level namespace
from .constructor import NORMS, LOCATORS, FORMATTERS, SCALES, PROJS
_globals = globals()
for _src in (NORMS, LOCATORS, FORMATTERS, SCALES, PROJS):
    for _key, _cls in _src.items():
        if isinstance(_cls, type):  # i.e. not a scale preset
            _globals[_cls.__name__] = _cls  # may overwrite proplot names

# Register objects
from .config import register_cmaps, register_cycles, register_colors, register_fonts
with benchmarks._benchmark('cmaps'):
    register_cmaps(default=True)
with benchmarks._benchmark('cycles'):
    register_cycles(default=True)
with benchmarks._benchmark('colors'):
    register_colors(default=True)
with benchmarks._benchmark('fonts'):
    register_fonts(default=True)

# Validate colormap names and propagate 'cycle' to 'axes.prop_cycle'
# NOTE: cmap.sequential also updates siblings 'cmap' and 'image.cmap'
from .config import rc
rcsetup.VALIDATE_REGISTERED_CMAPS = True
for _key in ('cycle', 'cmap.sequential', 'cmap.diverging', 'cmap.cyclic', 'cmap.qualitative'):  # noqa: E501
    try:
        rc[_key] = rc[_key]
    except ValueError as err:
        warnings._warn_proplot(f'Invalid user rc file setting: {err}')
        rc[_key] = 'Greys'  # fill value

# Validate color names now that colors are registered
from .config import rc_proplot, rc_matplotlib
rcsetup.VALIDATE_REGISTERED_COLORS = True
for _src in (rc_proplot, rc_matplotlib):
    for _key in _src:  # loop through unsynced properties
        if 'color' not in _key:
            continue
        # Likely has a color validator or derivative thereof; if not, harmless
        try:
            _src[_key] = _src[_key]
        except ValueError as err:
            warnings._warn_proplot(f'Invalid user rc file setting: {err}')
            _src[_key] = 'black'  # fill value
