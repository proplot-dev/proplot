#!/usr/bin/env python3
"""
Utilities for modifying proplot docstrings.
"""
import inspect
import re

import matplotlib.axes as maxes
import matplotlib.figure as mfigure
from matplotlib import rcParams as rc_matplotlib


class _SnippetManager(dict):
    """
    A simple database for snippets.
    """
    def __call__(self, obj):
        """
        Add snippets to the string or object using ``%(name)s`` substitution.
        This lets snippet keys have invalid variable names.
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
        Populate input strings with other snippets. Developers should take
        care to import modules in the correct order.
        """
        value = self(value)
        value = value.strip('\n')
        super().__setitem__(key, value)


def _obfuscate_signature(func):
    """
    Obfuscate a misleading or confusing call signature. Instead users
    should inspect the parameter table.
    """
    # Obfuscate signature by converting to *args **kwargs. Note this does
    # not change behavior of function! Copy parameters from a dummy function
    # because I'm too lazy to figure out inspect.Parameters API
    # See: https://stackoverflow.com/a/33112180/4970632
    sig = inspect.signature(func)
    sig_obfuscated = inspect.signature(lambda *args, **kwargs: None)
    func.__signature__ = sig.replace(
        parameters=tuple(sig_obfuscated.parameters.values())
    )
    return func


def _concatenate_original(func, prepend_summary=False):
    """
    Concatenate docstrings from a matplotlib axes method with a ProPlot axes
    method and obfuscate the call signature.
    """
    # NOTE: Do not bother inheriting from cartopy GeoAxes. Cartopy completely
    # truncates the matplotlib docstrings (which is kind of not great).
    # Get matplotlib axes func
    qual = func.__qualname__
    if 'Axes' in qual:
        cls = maxes.Axes
    elif 'Figure' in qual:
        cls = mfigure.Figure
    else:
        raise ValueError(f'Unexpected method {qual!r}. Must be Axes or Figure method.')
    doc = inspect.getdoc(func) or ''  # also dedents
    func_orig = getattr(cls, func.__name__, None)
    doc_orig = inspect.getdoc(func_orig)
    if not doc_orig:  # should never happen
        return func

    # Prepend summary
    if prepend_summary:
        regex = re.search(r'\.( | *\n|\Z)', doc_orig)
        if regex:
            doc = doc_orig[:regex.start() + 1] + '\n\n' + doc

    # Exit if this is documentation generated for website
    if rc_matplotlib['docstring.hardcopy']:
        func.__doc__ = doc
        return func

    # Concatenate docstrings and copy summary
    # Make sure different sections are very visible
    doc = f"""
=====================
ProPlot documentation
=====================

{doc}

========================
Matplotlib documentation
========================

{doc_orig}
"""
    func.__doc__ = inspect.cleandoc(doc)  # dedents and trims whitespace
    func = _obfuscate_signature(func)

    # Return
    return func


# Initiate snippets database
_snippet_manager = _SnippetManager()
