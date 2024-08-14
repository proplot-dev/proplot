#!/usr/bin/env python3
"""
Utilities for internal warnings and deprecations.
"""
import functools
import re
import sys
import warnings

from . import ic  # noqa: F401

# Internal modules omitted from warning message
REGEX_INTERNAL = re.compile(r"\A(matplotlib|mpl_toolkits|proplot)\.")

# Trivial warning class meant only to communicate the source of the warning
ProplotWarning = type("ProplotWarning", (UserWarning,), {})

# Add due to overwriting the module name
catch_warnings = warnings.catch_warnings
simplefilter = warnings.simplefilter


def next_release():
    """
    message indicating the next major release.
    """
    from .. import __version__

    try:
        # Find the first digit in the version string
        version_start = next(i for i, c in enumerate(__version__) if c.isdigit())
        num = int(__version__[version_start]) + 1
    except (StopIteration, ValueError, TypeError):
        string = "the next major release"
    else:
        which = "first" if num == 1 else "next"
        string = f"the {which} major release (version {num}.0.0)"
    return string


def _warn_proplot(message):
    """
    Emit a `ProplotWarning` and show the stack level outside of matplotlib and
    proplot. This is adapted from matplotlib's warning system.
    """
    frame = sys._getframe()
    stacklevel = 1
    while frame is not None:
        if not REGEX_INTERNAL.match(frame.f_globals.get("__name__", "")):
            break  # this is the first external frame
        frame = frame.f_back
        stacklevel += 1
    warnings.warn(message, ProplotWarning, stacklevel=stacklevel)


def _rename_objs(version, **kwargs):
    """
    Emit a basic deprecation warning after renaming function(s), method(s), or
    class(es). Each key should be an old name, and each argument should be the new
    object to point to. Do not document the deprecated object(s) to discourage use.
    """
    objs = []
    for old_name, new_obj in kwargs.items():
        new_name = new_obj.__name__
        message = (
            f"{old_name!r} was deprecated in version {version} and may be "
            f"removed in {next_release()}. Please use {new_name!r} instead."
        )
        if isinstance(new_obj, type):

            class _deprecated_class(new_obj):
                def __init__(self, *args, new_obj=new_obj, message=message, **kwargs):
                    _warn_proplot(message)
                    super().__init__(*args, **kwargs)

            _deprecated_class.__name__ = old_name
            objs.append(_deprecated_class)
        elif callable(new_obj):

            def _deprecated_function(*args, new_obj=new_obj, message=message, **kwargs):
                _warn_proplot(message)
                return new_obj(*args, **kwargs)

            _deprecated_function.__name__ = old_name
            objs.append(_deprecated_function)
        else:
            raise ValueError(f"Invalid deprecated object replacement {new_obj!r}.")
    if len(objs) == 1:
        return objs[0]
    else:
        return tuple(objs)


def _rename_kwargs(version, **kwargs_rename):
    """
    Emit a basic deprecation warning after removing or renaming keyword argument(s).
    Each key should be an old keyword, and each argument should be the new keyword
    or *instructions* for what to use instead.
    """

    def _decorator(func_orig):
        @functools.wraps(func_orig)
        def _deprecate_kwargs_wrapper(*args, **kwargs):
            for key_old, key_new in kwargs_rename.items():
                if key_old not in kwargs:
                    continue
                value = kwargs.pop(key_old)
                if key_new.isidentifier():
                    # Rename argument
                    kwargs[key_new] = value
                elif "{}" in key_new:
                    # Nice warning message, but user's desired behavior fails
                    key_new = key_new.format(value)
                _warn_proplot(
                    f"Keyword {key_old!r} was deprecated in version {version} and may "
                    f"be removed in {next_release()}. Please use {key_new!r} instead."
                )
            return func_orig(*args, **kwargs)

        return _deprecate_kwargs_wrapper

    return _decorator
