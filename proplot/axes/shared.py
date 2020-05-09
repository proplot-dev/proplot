"""
Shared utilities for axes classes.
"""
import functools


def _disable_decorator(msg):
    """
    Return a decorator that disables methods with message `msg`. The
    docstring is set to ``None`` so the ProPlot fork of automodapi doesn't add
    these methods to the website documentation. Users can still call
    help(ax.method) because python looks for superclass method docstrings if a
    docstring is empty.
    """
    def decorator(func):
        @functools.wraps(func)
        def _wrapper(self, *args, **kwargs):
            raise RuntimeError(msg.format(func.__name__))
        _wrapper.__doc__ = None
        return _wrapper
    return decorator
