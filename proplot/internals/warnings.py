#!/usr/bin/env python3
"""
Custom warning style and deprecation functions.
"""
import warnings
import functools


def _format_warning(message, category, filename, lineno, line=None):  # noqa: U100, E501
    """
    Warnings format monkey patch for warnings issued by ProPlot. See the
    `internal warning call signature \
<https://docs.python.org/3/library/warnings.html#warnings.showwarning>`__
    and the `default warning source code \
<https://github.com/python/cpython/blob/master/Lib/warnings.py>`__.
"""
    return f'{filename}:{lineno}: ProPlotWarning: {message}\n'  # needs newline


def _warn_proplot(message):
    """
    Temporarily apply the `_format_warning` monkey patch and emit the
    warning. Do not want to affect warnings emitted by other modules.
    """
    from . import _set_state
    with _set_state(warnings, formatwarning=_format_warning):
        warnings.warn(message, stacklevel=2)


def _rename_obj(old_name, new_obj, version=None):
    """
    Emit a basic deprecation warning after removing or renaming a function,
    method, or class. Do not document the deprecated object to discourage use.
    """
    new_name = new_obj.__name__
    version = 'a future version' if version is None else f'version {version}'

    def obj(*args, **kwargs):
        _warn_proplot(
            f'{old_name!r} is deprecated and will be removed in '
            f'{version}. Please use {new_name!r} instead.'
        )
        return new_obj(*args, **kwargs)
    obj.__name__ = old_name
    return obj


def _rename_kwargs(version=None, **kwargs_rename):
    """
    Emit a basic deprecation warning after removing or renaming function
    keyword arguments.
    """
    version = 'a future version' if version is None else f'version {version}'

    def decorator(func_orig):
        @functools.wraps(func_orig)
        def func(*args, **kwargs):
            for key_old, key_new in kwargs_rename.items():
                if key_old in kwargs:
                    if key_new is None:
                        del kwargs[key_old]
                        _warn_proplot(
                            f'Ignoring keyword arg {key_old!r}. This argument '
                            'is deprecated. Using it will raise an error '
                            'in {version}.'
                        )
                    else:
                        kwargs[key_new] = kwargs.pop(key_old)
                        _warn_proplot(
                            f'Keyword arg {key_old!r} is deprecated and will be '
                            f'removed in {version}. Please use {key_new!r} '
                            'instead.'
                        )
            return func_orig(*args, **kwargs)
        return func
    return decorator
