#!/usr/bin/env python3
"""
Utilities used internally by proplot.
"""
from . import defaults, docstring, timers, warnings  # noqa: F401


def _not_none(*args, default=None, **kwargs):
    """
    Return the first non-``None`` value. This is used with keyword arg
    aliases and for setting default values. Ugly name but clear purpose. Pass
    the `names` keyword arg to issue warning if multiple args were passed. Must
    be list of non-empty strings.
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


class _set_state(object):
    """
    Temporarily modify attribute(s) for an arbitrary object.
    """
    def __init__(self, obj, **kwargs):
        self._obj = obj
        self._kwargs = kwargs
        self._kwargs_orig = {
            key: getattr(obj, key) for key in kwargs if hasattr(obj, key)
        }

    def __enter__(self):
        for key, value in self._kwargs.items():
            setattr(self._obj, key, value)

    def __exit__(self, *args):  # noqa: U100
        for key in self._kwargs.keys():
            if key in self._kwargs_orig:
                setattr(self._obj, key, self._kwargs_orig[key])
            else:
                delattr(self._obj, key)
