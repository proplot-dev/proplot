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
    type_ = 'class' if isinstance(new_obj, type) else 'function'

    def obj(*args, **kwargs):
        _warn_proplot(
            f'{old_name!r} is deprecated and will be removed in {version}. '
            f'Please use {new_name!r} instead.',
        )
        return new_obj(*args, **kwargs)  # call function or instantiate class
    obj.__name__ = old_name
    obj.__doc__ = f"""
{type_.title()} {old_name!r} is deprecated and will be removed in {version}.
Please use {new_name!r} instead.
"""
    return obj


def _rename_kwargs(version=None, ignore=False, **kwargs_rename):
    """
    Emit a basic deprecation warning after removing or renaming function
    keyword arguments. Each key should be an old keyword, and each arguments
    should be the new keyword or a tuple of new keyword options.
    """
    version = 'a future version' if version is None else f'version {version}'

    def decorator(func_orig):
        @functools.wraps(func_orig)
        def func(*args, **kwargs):
            for key_old, key_new in kwargs_rename.items():
                if key_old in kwargs:
                    if ignore or not isinstance(key_new, str):
                        del kwargs[key_old]
                        message = f'Ignoring deprecated keyword arg {key_old!r}.'
                    else:
                        kwargs[key_new] = kwargs.pop(key_old)
                        message = (
                            f'Keyword arg {key_old!r} is deprecated and will '
                            f'be removed in {version}.'
                        )
                    if isinstance(key_new, str):
                        alternative = repr(key_new)
                    else:
                        alternative = ', '.join(map(repr, key_new))
                    _warn_proplot(f'{message} Please use {alternative} instead.')
            return func_orig(*args, **kwargs)
        return func
    return decorator
