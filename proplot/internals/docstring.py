#!/usr/bin/env python3
"""
Utilities for modifying proplot docstrings.
"""
import inspect
import re

import matplotlib.axes as maxes
import matplotlib.figure as mfigure
from matplotlib import rcParams as rc_matplotlib

#: A dictionary of docstring snippets.
snippets = {}


def add_snippets(obj):
    """
    A decorator that dedents docstrings with `inspect.getdoc` and adds
    un-indented snippets from the global `snippets` dictionary. This function
    uses ``%(name)s`` substitution rather than `str.format` substitution so
    that the `snippets` keys don't have to be valid variable names.
    """
    parts = {key: value.strip() for key, value in snippets.items()}
    if isinstance(obj, str):
        obj %= parts  # add snippets to a string
    else:
        obj.__doc__ = inspect.getdoc(obj)  # also dedents the docstring
        if obj.__doc__:
            obj.__doc__ %= parts  # insert snippets after dedent
    return obj


def obfuscate_signature(func):
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


def concatenate_original(func, prepend_summary=False):
    """
    Concatenate docstrings from a matplotlib axes method with a ProPlot axes
    method and obfuscate the call signature.
    """
    # NOTE: Originally had idea to use numpydoc.docscrape.NumpyDocString to
    # interpolate docstrings but *enormous* number of assupmtions would go into
    # this. And simple is better than complex.
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
    func = obfuscate_signature(func)

    # Return
    return func
