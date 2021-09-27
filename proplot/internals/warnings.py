#!/usr/bin/env python3
"""
Utilities for internal warnings and deprecations.
"""
import functools
import re
import sys
import warnings

# Internal modules omitted from warning message
REGEX_INTERNAL = re.compile(r'\A(matplotlib|mpl_toolkits|proplot)\.')

# Trivial warning class meant only to communicate the source of the warning
ProPlotWarning = type('ProPlotWarning', (UserWarning,), {})

# Add due to overwriting the module name
catch_warnings = warnings.catch_warnings
simplefilter = warnings.simplefilter


def _warn_proplot(message, action=None):
    """
    Emit a `ProPlotWarning` and show the stack level outside of matplotlib and
    proplot. This is adapted from matplotlib's warning system.
    """
    frame = sys._getframe()
    stacklevel = 1
    while frame is not None:
        if not REGEX_INTERNAL.match(frame.f_globals.get('__name__', '')):
            break  # this is the first external frame
        frame = frame.f_back
        stacklevel += 1
    with warnings.catch_warnings():
        if action:  # used internally
            warnings.simplefilter(action, ProPlotWarning)
        warnings.warn(message, ProPlotWarning, stacklevel=stacklevel)


def _rename_objs(version, **kwargs):
    """
    Emit a basic deprecation warning after renaming function(s), method(s), or
    class(es). Each key should be an old name, and each argument should be the new
    object to point to. Do not document the deprecated object(s) to discourage use.
    """
    wrappers = []
    for old_name, new_obj in kwargs.items():
        new_name = new_obj.__name__
        message = (
            f'{old_name!r} was deprecated in version {version} and will be '
            f'removed in a future release. Please use {new_name!r} instead.'
        )
        if isinstance(new_obj, type):
            class _deprecate_obj(new_obj):
                def __init__(self, *args, new_obj=new_obj, message=message, **kwargs):
                    _warn_proplot(message)
                    super().__init__(*args, **kwargs)
        elif callable(new_obj):
            def _deprecate_obj(*args, new_obj=new_obj, message=message, **kwargs):
                _warn_proplot(message)
                return new_obj(*args, **kwargs)
        else:
            raise ValueError(f'Invalid deprecated object replacement {new_obj!r}.')
        _deprecate_obj.__name__ = old_name
        wrappers.append(_deprecate_obj)
    if len(wrappers) == 1:
        return wrappers[0]
    else:
        return tuple(wrappers)


def _rename_kwargs(version, **kwargs_rename):
    """
    Emit a basic deprecation warning after removing or renaming keyword argument(s).
    Each key should be an old keyword, and each argument should be the new keyword
    or *instructions* for what to use instead.
    """
    def decorator(func_orig):
        @functools.wraps(func_orig)
        def _deprecate_kwargs(*args, **kwargs):
            for key_old, key_new in kwargs_rename.items():
                if key_old not in kwargs:
                    continue
                value = kwargs.pop(key_old)
                if key_new.isidentifier():
                    # Rename argument
                    kwargs[key_new] = value
                elif '{}' in key_new:
                    # Nice warning message, but user's desired behavior fails
                    key_new = key_new.format(value)
                _warn_proplot(
                    f'Keyword {key_old!r} was deprecated in version {version} and will '
                    f'be removed in a future release. Please use {key_new!r} instead.'
                )
            return func_orig(*args, **kwargs)
        return _deprecate_kwargs
    return decorator
