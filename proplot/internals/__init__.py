#!/usr/bin/env python3
"""
Utilities used internally by proplot.
"""
from . import rcsetup, docstring, timers, warnings  # noqa: F401
try:  # print debugging
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa


def _not_none(*args, default=None, **kwargs):
    """
    Return the first non-``None`` value. This is used with keyword arg aliases and
    for setting default values. Use `kwargs` to issue warnings when multiple
    non-``None`` values were passed. Use `args` to just ignore extra values.
    """
    first = default
    if args and kwargs:
        raise ValueError('_not_none can only be used with args or kwargs.')
    elif args:
        for arg in args:
            if arg is not None:
                first = arg
                break
    elif kwargs:
        for name, arg in list(kwargs.items()):
            if arg is not None:
                first = arg
                break
        kwargs = {name: arg for name, arg in kwargs.items() if arg is not None}
        if len(kwargs) > 1:
            warnings._warn_proplot(
                f'Got conflicting or duplicate keyword args: {kwargs}. '
                'Using the first one.'
            )
    return first


class _dummy_context(object):
    """
    Dummy context manager.
    """
    def __init__(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args):  # noqa: U100
        pass


class _set_state(object):
    """
    Temporarily modify attribute(s) for an arbitrary object.
    """
    def __init__(self, obj, **kwargs):
        self._obj = obj
        self._attrs_new = kwargs
        self._attrs_prev = {
            key: getattr(obj, key) for key in kwargs if hasattr(obj, key)
        }

    def __enter__(self):
        for key, value in self._attrs_new.items():
            setattr(self._obj, key, value)

    def __exit__(self, *args):  # noqa: U100
        for key in self._attrs_new.keys():
            if key in self._attrs_prev:
                setattr(self._obj, key, self._attrs_prev[key])
            else:
                delattr(self._obj, key)


class _version(list):
    """
    Casual version parser for MAJOR.MINOR version strings. Do not want to add
    'packaging' dependency and only care about major and minor tags.
    """
    def __repr__(self):
        return f'version({self._version})'

    def __init__(self, version):
        try:
            major, minor, *_ = version.split('.')
            major, minor = int(major), int(minor)
        except (ValueError, AttributeError):
            warnings._warn_proplot(
                f"Invalid version {version!r}. Defaulting to '0.0'."
            )
            major = minor = 0
        self._version = version
        super().__init__([major, minor])  # then use builtin python list sorting


# Add matplotlib and cartopy versions
import matplotlib as _
_version_mpl = _version(_.__version__)
try:
    import cartopy as _
    _version_cartopy = _version(_.__version__)
except ImportError:
    _version_cartopy = [0, 0]
