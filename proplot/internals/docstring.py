#!/usr/bin/env python3
"""
Utilities for modifying proplot docstrings.
"""
# WARNING: To check every docstring in the package for
# unfilled snippets simply use the following code:
# >>> import proplot as pplt
# ... seen = set()
# ... def _iter_doc(objs):
# ...     if objs in seen:
# ...         return
# ...     seen.add(objs)
# ...     for attr in dir(objs):
# ...         obj = getattr(objs, attr, None)
# ...         if callable(obj) and hasattr(obj, '__doc__'):
# ...             if obj in seen:
# ...                 continue
# ...             seen.add(obj)
# ...             if obj.__doc__ and '%(' in obj.__doc__:
# ...                 yield obj.__name__
# ...             yield from _iter_doc(obj)
# ... print(*_iter_doc(pplt))
import inspect
import re

import matplotlib.axes as maxes
import matplotlib.figure as mfigure
from matplotlib import rcParams as rc_matplotlib

from . import ic  # noqa: F401


def _obfuscate_kwargs(func):
    """
    Obfuscate keyword args.
    """
    return _obfuscate_signature(func, lambda **kwargs: None)


def _obfuscate_params(func):
    """
    Obfuscate all parameters.
    """
    return _obfuscate_signature(func, lambda *args, **kwargs: None)


def _obfuscate_signature(func, dummy):
    """
    Obfuscate a misleading or incomplete call signature.
    Instead users should inspect the parameter table.
    """
    # Obfuscate signature by converting to *args **kwargs. Note this does
    # not change behavior of function! Copy parameters from a dummy function
    # because I'm too lazy to figure out inspect.Parameters API
    # See: https://stackoverflow.com/a/33112180/4970632
    sig = inspect.signature(func)
    sig_repl = inspect.signature(dummy)
    func.__signature__ = sig.replace(parameters=tuple(sig_repl.parameters.values()))
    return func


def _concatenate_inherited(func, prepend_summary=False):
    """
    Concatenate docstrings from a matplotlib axes method with a proplot
    axes method and obfuscate the call signature.
    """
    # Get matplotlib axes func
    # NOTE: Do not bother inheriting from cartopy GeoAxes. Cartopy completely
    # truncates the matplotlib docstrings (which is kind of not great).
    qual = func.__qualname__
    if "Axes" in qual:
        cls = maxes.Axes
    elif "Figure" in qual:
        cls = mfigure.Figure
    else:
        raise ValueError(f"Unexpected method {qual!r}. Must be Axes or Figure method.")
    doc = inspect.getdoc(func) or ""  # also dedents
    func_orig = getattr(cls, func.__name__, None)
    doc_orig = inspect.getdoc(func_orig)
    if not doc_orig:  # should never happen
        return func

    # Optionally prepend the function summary
    # Concatenate docstrings only if this is not generated for website
    regex = re.search(r"\.( | *\n|\Z)", doc_orig)
    if regex and prepend_summary:
        doc = doc_orig[: regex.start() + 1] + "\n\n" + doc
    if not rc_matplotlib["docstring.hardcopy"]:
        doc = f"""
=====================
Proplot documentation
=====================

{doc}

========================
Matplotlib documentation
========================

{doc_orig}
"""

    # Return docstring
    # NOTE: Also obfuscate parameters to avoid partial coverage of call signatures
    func.__doc__ = inspect.cleandoc(doc)
    func = _obfuscate_params(func)
    return func


class _SnippetManager(dict):
    """
    A simple database for handling documentation snippets.
    """

    def __call__(self, obj):
        """
        Add snippets to the string or object using ``%(name)s`` substitution. Here
        ``%(name)s`` is used rather than ``.format`` to support invalid identifiers.
        """
        if isinstance(obj, str):
            obj %= self  # add snippets to a string
        else:
            obj.__doc__ = inspect.getdoc(obj)  # also dedents the docstring
            if obj.__doc__:
                obj.__doc__ %= self  # insert snippets after dedent
        return obj

    def __setitem__(self, key, value):
        """
        Populate input strings with other snippets and strip newlines. Developers
        should take care to import modules in the correct order.
        """
        value = self(value)
        value = value.strip("\n")
        super().__setitem__(key, value)


# Initiate snippets database
_snippet_manager = _SnippetManager()
