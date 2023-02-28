#!/usr/bin/env python3
"""
Utilities for global styles and stylesheets.
"""
import os

import matplotlib as mpl
import matplotlib.style as mstyle
from matplotlib import RcParams
from matplotlib import rcParamsDefault as _rc_matplotlib_default
from matplotlib import rcParamsOrig as _rc_matplotlib_original

from . import ic  # noqa: F401
from . import defaults, warnings

# Matplotlib stylesheets
# NOTE: The 'proplot' style is not registered in matplotlib but instead is specially
# handled by proplot.config.use_style (similar to 'default' in matplotlib.style.use)
STYLE_ALIASES = {
    None: 'proplot',
    '538': 'fivethirtyeight',
    'mpl15': 'classic',
    'mpl20': 'default',
    'matplotlib': 'default',
}


def _filter_style_dict(kw, warn_blacklisted=True):
    """
    Filter out blacklisted style parameters.
    """
    # NOTE: This implements bugfix https://github.com/matplotlib/matplotlib/pull/17252
    # critical for proplot because we always run style.use() when the configurator is
    # initialized. Without fix backend resets every time you import proplot. Otherwise
    # this is just a copy of _remove_blacklisted_params in matplotlib.style.core.
    kw_filtered = {}
    for key in kw:
        if key in mstyle.core.STYLE_BLACKLIST:
            if warn_blacklisted:
                warnings._warn_proplot(
                    f'Dictionary includes a parameter, {key!r}, that is '
                    'not related to style. Ignoring.'
                )
        else:
            kw_filtered[key] = kw[key]
    return kw_filtered


def _get_style_dict(style, **kwargs):
    """
    Get the default or original rc dictionary with deprecated parameters filtered.
    """
    # WARNING: Some deprecated rc params remain in dictionary as None in matplotlib
    # < 3.5 (?) so manually filter them out (_deprecated_set is matplotlib < 3.4 and
    # examples.directory was handled specially inside RcParams in matplotlib < 3.2).
    style = STYLE_ALIASES.get(style, style)
    extras = ('default', 'original', 'proplot', 'texgyre', *STYLE_ALIASES)
    if style in ('default', 'original', 'proplot'):
        kw = _rc_matplotlib_original if style == 'original' else _rc_matplotlib_default
        kw = _filter_style_dict(kw, warn_blacklisted=False)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', mpl.MatplotlibDeprecationWarning)
            kw = dict(RcParams(kw))  # filter and translate deprecations
        for attr in ('_deprecated_set', '_deprecated_remain_as_none'):
            deprecated = getattr(mpl, attr, ())
            for key in deprecated:
                kw.pop(key, None)  # remove
        kw.pop('examples.directory', None)  # special case for matplotlib < 3.2
    elif style == 'texgyre':
        kw = {
            'font.' + family: defaults._rc_matplotlib_override['font.' + family]
            for family in ('serif', 'sans-serif', 'monospace', 'cursive', 'fantasy')
        }
    elif style in mstyle.library:
        kw = mstyle.library[style]
    else:
        try:
            kw = mpl.rc_params_from_file(style, use_default_template=False)
        except IOError:
            raise IOError(
                f'Style {style!r} not found in the style library and input '
                'is not a valid URL or file path. Available styles are: '
                + ', '.join(map(repr, (*extras, *mstyle.library)))
                + '.'
            )
    return _filter_style_dict(kw, **kwargs)


def _infer_style_dict(kw):
    """
    Infer proplot parameter values from matploglib stylesheet parameters.
    """
    # NOTE: This helps make native matplotlib stylesheets extensible to proplot
    # figures without crazy font size and color differences.
    kw_proplot = {}
    mpl_to_proplot = {
        'xtick.labelsize': (
            'tick.labelsize', 'grid.labelsize',
        ),
        'ytick.labelsize': (
            'tick.labelsize', 'grid.labelsize',
        ),
        'axes.titlesize': (
            'abc.size', 'suptitle.size', 'title.size',
            'leftlabel.size', 'rightlabel.size',
            'toplabel.size', 'bottomlabel.size',
        ),
        'text.color': (
            'abc.color', 'suptitle.color', 'title.color',
            'tick.labelcolor', 'grid.labelcolor',
            'leftlabel.color', 'rightlabel.color',
            'toplabel.color', 'bottomlabel.color',
        ),
    }
    for key, params in mpl_to_proplot.items():
        if key in kw:
            value = kw[key]
            for param in params:
                kw_proplot[param] = value
    return kw_proplot


def _parse_style_spec(styles, allow_dictionary=True, **kwargs):
    """
    Validate the matplotlib style or style list.
    """
    # NOTE: Hold off on filtering out non-style parameters. Will issue warnings
    # on subsequent assignment inside 'use_style'.
    # NOTE: Unlike matplotlib this begins all styles from the custom base state of
    # the matplotlib settings + proplot fonts.
    kw_matplotlib = {}
    if styles is None or isinstance(styles, (str, dict, os.PathLike)):
        styles = [styles]
    else:
        styles = list(styles)
    if 'classic' in styles:
        styles = ['default', *styles]
    elif any(s in styles for s in ('default', 'original', 'proplot')):
        pass
    else:
        styles = ['default', 'texgyre', *styles]
    for style in styles:
        if allow_dictionary and isinstance(style, dict):
            kw = style
        else:
            kw = _get_style_dict(style, **kwargs)
        kw_matplotlib.update(kw)
    kw_proplot = _infer_style_dict(kw_matplotlib)
    kw_proplot['style'] = tuple(styles)
    return kw_matplotlib, kw_proplot
