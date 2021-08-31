#!/usr/bin/env python3
"""
The starting point for creating ProPlot figures.
"""
import matplotlib.pyplot as plt

from . import figure as pfigure
from . import gridspec as pgridspec
from .internals import ic  # noqa: F401
from .internals import _not_none, _pop_params, docstring

__all__ = [
    'figure',
    'subplots',
    'show',
    'close',
    'switch_backend',
    'ion',
    'ioff',
    'isinteractive',
]


# Docstrings
_pyplot_docstring = """
This is included so you don't have to import `~matplotlib.pyplot`.
"""
docstring._snippet_manager['ui.pyplot'] = _pyplot_docstring


def _parse_figsize(kwargs):
    """
    Translate `figsize` into proplot-specific `figwidth` and `figheight` keys.
    """
    # WARNING: Cannot have Figure.__init__() interpret figsize() because
    # the figure manager fills it with the matplotlib default.
    figsize = kwargs.pop('figsize', None)
    figwidth = kwargs.pop('figwidth', None)
    figheight = kwargs.pop('figheight', None)
    if figsize is not None:
        figsize_width, figsize_height = figsize
        figwidth = _not_none(figwidth=figwidth, figsize_width=figsize_width)
        figheight = _not_none(figheight=figheight, figsize_height=figsize_height)
    kwargs['figwidth'] = figwidth
    kwargs['figheight'] = figheight


@docstring._snippet_manager
def show(*args, **kwargs):
    """
    Call `matplotlib.pyplot.show`.
    %(ui.pyplot)s

    Parameters
    ----------
    *args, **kwargs
        Passed to `matplotlib.pyplot.show`.
    """
    return plt.show(*args, **kwargs)


@docstring._snippet_manager
def close(*args, **kwargs):
    """
    Call `matplotlib.pyplot.close`.
    %(ui.pyplot)s

    Parameters
    ----------
    *args, **kwargs
        Passed to `matplotlib.pyplot.close`.
    """
    return plt.close(*args, **kwargs)


@docstring._snippet_manager
def switch_backend(*args, **kwargs):
    """
    Call `matplotlib.pyplot.switch_backend`.
    %(ui.pyplot)s

    Parameters
    ----------
    *args, **kwargs
        Passed to `matplotlib.pyplot.switch_backend`.
    """
    return plt.switch_backend(*args, **kwargs)


@docstring._snippet_manager
def ion():
    """
    Call `matplotlib.pyplot.ion`.
    %(ui.pyplot)s
    """
    return plt.ion()


@docstring._snippet_manager
def ioff():
    """
    Call `matplotlib.pyplot.ioff`.
    %(ui.pyplot)s
    """
    return plt.ioff()


@docstring._snippet_manager
def isinteractive():
    """
    Call `matplotlib.pyplot.isinteractive`.
    %(ui.pyplot)s
    """
    return plt.isinteractive()


@docstring._snippet_manager
def figure(**kwargs):
    """
    Create an empty figure. Subplots can be subsequently added using
    `~proplot.figure.Figure.add_subplot` or `~proplot.figure.Figure.subplots`.
    This command is analogous to `matplotlib.pyplot.figure`.

    Parameters
    ----------
    %(figure.figure)s

    Other parameters
    ----------------
    **kwargs
        Passed to `matplotlib.figure.Figure`.

    See also
    --------
    proplot.ui.subplots
    proplot.figure.Figure.add_subplot
    proplot.figure.Figure.subplots
    proplot.figure.Figure
    matplotlib.figure.Figure
    """
    _parse_figsize(kwargs)
    return plt.figure(FigureClass=pfigure.Figure, **kwargs)


@docstring._snippet_manager
def subplots(*args, **kwargs):
    """
    Create a figure with a single subplot or an arbitrary grid of
    subplots. Uses `figure` and `proplot.figure.Figure.subplots`.
    This command is analogous to `matplotlib.pyplot.subplots`.

    Parameters
    ----------
    %(figure.subplots_params)s

    Other parameters
    ----------------
    %(figure.figure)s
    **kwargs
        Passed to `matplotlib.figure.Figure`.

    Returns
    -------
    fig : `proplot.figure.Figure`
        The figure instance.
    axs : `proplot.gridspec.SubplotGrid`
        The axes instances stored in a `~proplot.gridspec.SubplotGrid`.

    See also
    --------
    proplot.ui.figure
    proplot.figure.Figure.subplots
    proplot.gridspec.SubplotGrid
    proplot.figure.Figure
    matplotlib.figure.Figure
    """
    _parse_figsize(kwargs)
    kw_subplots = _pop_params(kwargs, pfigure.Figure.add_subplots)
    kw_gridspec = _pop_params(kwargs, pgridspec.GridSpec._update_params)
    fig = figure(**kwargs)
    axs = fig.add_subplots(*args, **kw_subplots, **kw_gridspec)
    return fig, axs
