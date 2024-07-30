#!/usr/bin/env python3
"""
The starting point for creating proplot figures.
"""
import matplotlib.pyplot as plt

from . import axes as paxes
from . import figure as pfigure
from . import gridspec as pgridspec
from .internals import ic  # noqa: F401
from .internals import _not_none, _pop_params, _pop_props, _pop_rc, docstring

__all__ = [
    "figure",
    "subplot",
    "subplots",
    "show",
    "close",
    "switch_backend",
    "ion",
    "ioff",
    "isinteractive",
]


# Docstrings
_pyplot_docstring = """
This is included so you don't have to import `~matplotlib.pyplot`.
"""
docstring._snippet_manager["ui.pyplot"] = _pyplot_docstring


def _parse_figsize(kwargs):
    """
    Translate `figsize` into proplot-specific `figwidth` and `figheight` keys.
    """
    # WARNING: Cannot have Figure.__init__() interpret figsize() because
    # the figure manager fills it with the matplotlib default.
    figsize = kwargs.pop("figsize", None)
    figwidth = kwargs.pop("figwidth", None)
    figheight = kwargs.pop("figheight", None)
    if figsize is not None:
        figsize_width, figsize_height = figsize
        figwidth = _not_none(figwidth=figwidth, figsize_width=figsize_width)
        figheight = _not_none(figheight=figheight, figsize_height=figsize_height)
    kwargs["figwidth"] = figwidth
    kwargs["figheight"] = figheight


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
        Passed to `proplot.figure.Figure.format`.

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
def subplot(**kwargs):
    """
    Return a figure and a single subplot.
    This command is analogous to `matplotlib.pyplot.subplot`,
    except the figure instance is also returned.

    Other parameters
    ----------------
    %(figure.figure)s
    **kwargs
        Passed to `proplot.figure.Figure.format` or the
        projection-specific ``format`` command for the axes.

    Returns
    -------
    fig : `proplot.figure.Figure`
        The figure instance.
    ax : `proplot.axes.Axes`
        The axes instance.

    See also
    --------
    proplot.ui.figure
    proplot.figure.Figure.subplot
    proplot.figure.Figure
    matplotlib.figure.Figure
    """
    _parse_figsize(kwargs)
    rc_kw, rc_mode = _pop_rc(kwargs)
    kwsub = _pop_props(kwargs, "patch")  # e.g. 'color'
    kwsub.update(_pop_params(kwargs, pfigure.Figure._parse_proj))
    for sig in paxes.Axes._format_signatures.values():
        kwsub.update(_pop_params(kwargs, sig))
    kwargs["aspect"] = kwsub.pop("aspect", None)  # keyword conflict
    fig = figure(rc_kw=rc_kw, **kwargs)
    ax = fig.add_subplot(rc_kw=rc_kw, **kwsub)
    return fig, ax


@docstring._snippet_manager
def subplots(*args, **kwargs):
    """
    Return a figure and an arbitrary grid of subplots.
    This command is analogous to `matplotlib.pyplot.subplots`,
    except the subplots are stored in a `~proplot.gridspec.SubplotGrid`.

    Parameters
    ----------
    %(figure.subplots_params)s

    Other parameters
    ----------------
    %(figure.figure)s
    **kwargs
        Passed to `proplot.figure.Figure.format` or the
        projection-specific ``format`` command for each axes.

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
    rc_kw, rc_mode = _pop_rc(kwargs)
    kwsubs = _pop_props(kwargs, "patch")  # e.g. 'color'
    kwsubs.update(_pop_params(kwargs, pfigure.Figure._add_subplots))
    kwsubs.update(_pop_params(kwargs, pgridspec.GridSpec._update_params))
    for sig in paxes.Axes._format_signatures.values():
        kwsubs.update(_pop_params(kwargs, sig))
    for key in ("subplot_kw", "gridspec_kw"):  # deprecated args
        if key in kwargs:
            kwsubs[key] = kwargs.pop(key)
    kwargs["aspect"] = kwsubs.pop("aspect", None)  # keyword conflict
    fig = figure(rc_kw=rc_kw, **kwargs)
    axs = fig.add_subplots(*args, rc_kw=rc_kw, **kwsubs)
    return fig, axs
