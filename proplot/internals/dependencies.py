#!/usr/bin/env python3
"""
Utilities for handling dependencies and version changes.
"""
from numbers import Real

from . import warnings


class _version(list):
    """
    Casual parser for ``major.minor`` style version strings. We do not want to
    add a 'packaging' dependency and only care about major and minor tags.
    """
    def __str__(self):
        return self._version

    def __repr__(self):
        return f'version({self._version})'

    def __init__(self, version):
        try:
            if isinstance(version, Real):
                version = str(version)
            if not isinstance(version, str):
                version = '.'.join(map(str, version))
            major, minor, *_ = version.split('.')
            major, minor = int(major or 0), int(minor or 0)
        except Exception:
            warnings._warn_proplot(
                f'Unexpected version {version!r}. Interpreting as 0.0.'
            )
            major = minor = 0
        self._version = version
        super().__init__((major, minor))  # then use builtin python list sorting

    def __eq__(self, other):
        return super().__eq__(_version(other))

    def __ne__(self, other):
        return super().__ne__(_version(other))

    def __gt__(self, other):
        return super().__gt__(_version(other))

    def __lt__(self, other):
        return super().__lt__(_version(other))

    def __ge__(self, other):
        return super().__ge__(_version(other))

    def __le__(self, other):
        return super().__le__(_version(other))


# Matplotlib version
import matplotlib  # isort:skip
_version_mpl = _version(matplotlib.__version__)

# Cartopy version
try:
    import cartopy
except ImportError:
    _version_cartopy = _version((0, 0))
else:
    _version_cartopy = _version(cartopy.__version__)
