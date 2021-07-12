#!/usr/bin/env python3
"""
Utilities for modifying proplot docstrings.
"""
import inspect

#: Dictionary of docstring snippets.
snippets = {}


def add_snippets(obj):
    """
    Decorator that dedents docstrings with `inspect.getdoc` and adds
    un-indented snippets from the global `snippets` dictionary. This function
    uses ``%(name)s`` substitution rather than `str.format` substitution so
    that the `snippets` keys can be invalid variable names.
    """
    parts = {key: value.strip() for key, value in snippets.items()}
    if isinstance(obj, str):
        obj %= parts
    else:
        obj.__doc__ = inspect.getdoc(obj)
        if obj.__doc__:
            obj.__doc__ %= parts
    return obj
