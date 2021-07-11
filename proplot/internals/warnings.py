#!/usr/bin/env python3
"""
Custom warning style and deprecation functions.
"""
import functools
import re
import sys
import warnings

ProPlotWarning = type('ProPlotWarning', (UserWarning,), {})


def _warn_proplot(message):
    """
    Emit a `ProPlotWarning` and show the stack level outside of matplotlib and proplot.
    """
    frame = sys._getframe()
    stacklevel = 1
    while True:
        if frame is None:
            break  # when called in embedded context may hit frame is None
        if not re.match(
            r'\A(matplotlib|mpl_toolkits|proplot)\.',
            frame.f_globals.get('__name__', '')
        ):
            break
        frame = frame.f_back
        stacklevel += 1
    warnings.warn(message, ProPlotWarning, stacklevel=stacklevel)


def _deprecate_getter_setter(version, property):
    """
    Generate `set_name` and `get_name` methods for property setters and getters,
    and issue warnings when they are used.
    """
    def getter(self):
        _warn_proplot(
            f'get_{property}() was deprecated in {version}. Please use '
            f'{type(self).__name__}.{property} instead.'
        )
        return getattr(self, '_' + property)

    def setter(self, value):
        _warn_proplot(
            f'set_{property}() was deprecated in {version}. The property is '
            f'now read-only.'
        )
        return

    getter.__name__ = f'get_{property}'
    setter.__name__ = f'set_{property}'

    return getter, setter


def _rename_objs(version, **kwargs):
    """
    Emit a basic deprecation warning after renaming function(s), method(s), or
    class(es). Do not document the deprecated object to discourage use.
    """
    wrappers = []
    for old_name, new_obj in kwargs.items():
        # Add info as keywords to avoid overwriting by loop scope
        def deprecate_obj(*args, _old_name=old_name, _new_obj=new_obj, **kwargs):
            _new_name = _new_obj.__name__
            _warn_proplot(
                f'{_old_name!r} was deprecated in version {version} and will be '
                f'removed in the next major release. Please use {_new_name!r} instead.'
            )
            return _new_obj(*args, **kwargs)
        # Replace name
        deprecate_obj.__name__ = old_name
        wrappers.append(deprecate_obj)
    if len(wrappers) == 1:
        return wrappers[0]
    else:
        return tuple(wrappers)


def _rename_kwargs(version, **kwargs_rename):
    """
    Emit a basic deprecation warning after removing or renaming function keyword
    arguments. Each key should be an old keyword, and each argument should be the
    new keyword or *instructions* for what to use instead.
    """
    def decorator(func_orig):
        @functools.wraps(func_orig)
        def deprecate_kwargs(*args, **kwargs):
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
                    f'Keyword arg {key_old!r} was deprecated in {version} and will be '
                    'removed in the next major release. Please use {key_new!r} instead.'
                )
            return func_orig(*args, **kwargs)
        return deprecate_kwargs
    return decorator
