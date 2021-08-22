#!/usr/bin/env python3
"""
Utilities used internally by proplot.
"""
import inspect

import numpy as np

from . import benchmarks, dependencies, docstring, rcsetup, warnings  # noqa: F401
from .dependencies import _version, _version_cartopy, _version_mpl  # noqa: F401
from .docstring import _snippet_manager  # noqa: F401

try:  # print debugging
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed
    ic = lambda *args: print(*args)  # noqa: E731

INTERNAL_PARAMS = {  # silently pop these if we don't reach certain internal utilities
    'to_centers',
    'line_plot',
    'contour_plot',
    'default_cmap',
    'default_discrete',
    'skip_autolev',
}

# Alias dictionaries. This package only works with a subset of available artists
# and keywords so we simply create our own system rather than working with
# matplotlib's normalize_kwargs and _alias_maps.
# WARNING: Add pseudo-props 'edgewidth' and 'fillcolor' for patch edges and faces
# WARNING: Critical that alias does not appear in key dict or else _translate_kwargs
# will overwrite settings with None after popping them!
_alias_dicts = {
    'rgba': {
        'red': ('r',),
        'green': ('g',),
        'blue': ('b',),
        'alpha': ('a',),
    },
    'hsla': {
        'hue': ('h',),
        'saturation': ('s', 'c', 'chroma'),
        'luminance': ('l',),
        'alpha': ('a',),
    },
    'line': {  # copied from lines.py but expanded to include plurals
        'antialiased': ('aa',),
        'alpha': ('a', 'alphas'),
        'color': ('c', 'colors'),
        'linewidth': ('lw', 'linewidths'),
        'linestyle': ('ls', 'linestyles'),
        'drawstyle': ('ds', 'drawstyles'),
        'dashes': ('d',),
        'marker': ('m', 'markers'),
        'markersize': ('s', 'ms', 'markersizes'),
        'markeredgecolor': ('mec', 'markeredgecolors'),
        'markeredgewidth': ('mew', 'markeredgewidths'),
        'markerfacecolor': ('mfc', 'markerfacecolors', 'mc', 'markercolor', 'markercolors'),  # noqa: E501
        'fillstyle': ('fs', 'fillstyles', 'mfs', 'markerfillstyle', 'markerfillstyles'),
        'zorder': ('z', 'zorders'),
    },
    'collection': {  # NOTE: face color is ignored for line collections
        'alphas': ('a', 'alpha'),
        'colors': ('c', 'color'),
        'edgecolors': ('ec', 'edgecolor'),
        'facecolors': ('fc', 'fillcolor', 'fillcolors'),
        'linewidths': ('lw', 'linewidth', 'ew', 'edgewidth', 'edgewidths'),
        'linestyles': ('ls', 'linestyle'),
        'zorder': ('z', 'zorders'),
    },
    'patch': {
        'alpha': ('a', 'alphas', 'facealpha', 'facealphas', 'fillalpha', 'fillalphas'),
        'color': ('c', 'colors'),
        'edgecolor': ('ec', 'edgecolors'),
        'facecolor': ('fc', 'facecolors', 'fillcolor', 'fillcolors'),
        'linewidth': ('lw', 'linewidths', 'ew', 'edgewidth', 'edgewidths'),
        'linestyle': ('ls', 'linestyles'),
        'zorder': ('z', 'zorders'),
        'hatch': ('h', 'hatching'),
    },
    'text': {
        'text': (),
        'color': ('c', 'fontcolor'),  # NOTE: see text.py source code
        'fontfamily': ('family',),
        'fontname': ('name',),
        'fontsize': ('size',),
        'fontstretch': ('stretch',),
        'fontstyle': ('style',),
        'fontvariant': ('variant',),
        'fontweight': ('weight',),
        'fontproperties': ('fp', 'font', 'font_properties'),
        'zorder': ('z', 'zorders'),
    },
}


def _not_none(*args, default=None, **kwargs):
    """
    Return the first non-``None`` value. This is used with keyword arg aliases and
    for setting default values. Use `kwargs` to issue warnings when multiple passed.
    """
    first = default
    if args and kwargs:
        raise ValueError('_not_none can only be used with args or kwargs.')
    elif args:
        for arg in args:
            if arg is not None:
                first = arg
                break
    elif kwargs:
        for name, arg in list(kwargs.items()):
            if arg is not None:
                first = arg
                break
        kwargs = {name: arg for name, arg in kwargs.items() if arg is not None}
        if len(kwargs) > 1:
            warnings._warn_proplot(
                f'Got conflicting or duplicate keyword args: {kwargs}. '
                'Using the first one.'
            )
    return first


def _keyword_to_positional(options, *args, allow_extra=False, **kwargs):
    """
    Translate keyword arguments to positional arguments. Permit omitted
    arguments so that plotting functions can infer values.
    """
    nargs, nopts = len(args), len(options)
    if nargs > nopts and not allow_extra:
        raise ValueError(f'Expected up to {nopts} positional arguments. Got {nargs}.')
    args = list(args)
    args.extend(None for _ in range(nopts - nargs))  # fill missing args
    for idx, keys in enumerate(options):
        if isinstance(keys, str):
            keys = (keys,)
        opts = {}
        if args[idx] is not None:  # positional args have first priority
            opts[keys[0] + '_positional'] = args[idx]
        for key in keys:  # keyword args
            opts[key] = kwargs.pop(key, None)
        args[idx] = _not_none(**opts)  # may reassign None
    return args, kwargs


def _translate_kwargs(input, output, *keys, **aliases):
    """
    The driver function.
    """
    aliases.update({key: () for key in keys})
    for key, aliases in aliases.items():
        aliases = (aliases,) if isinstance(aliases, str) else aliases
        opts = {key: input.pop(key, None) for key in (key, *aliases)}
        value = _not_none(**opts)
        if value is not None:
            output[key] = value
    return output


def _translate_props(input, output, *categories, prefix=None, ignore=None):  # noqa: E501
    """
    The driver function.
    """
    # Get properties
    prefix = prefix or ''  # e.g. 'box' for boxlw, boxlinewidth, etc.
    for category in categories:
        if category not in _alias_dicts:
            raise ValueError(f'Invalid alias category {category!r}.')
        for key, aliases in _alias_dicts[category].items():
            if isinstance(aliases, str):
                aliases = (aliases,)
            opts = {prefix + alias: input.pop(prefix + alias, None) for alias in (key, *aliases)}  # noqa: E501
            prop = _not_none(**opts)
            if prop is not None:
                output[key] = prop
    # Ignore properties (e.g., ignore 'marker' properties)
    ignore = ignore or ()
    if isinstance(ignore, str):
        ignore = (ignore,)
    for string in ignore:
        for key in tuple(output):
            if string in key:
                value = output.pop(key)
                warnings._warn_proplot(f'Ignoring property {key}={value!r}.')
    return output


def _pop_kwargs(src, *keys, **aliases):
    """
    Pop out input properties and return them in a new dictionary.
    """
    return _translate_kwargs(src, {}, *keys, **aliases)


def _process_kwargs(src, *keys, **aliases):
    """
    Translate input properties and add translated names to the original dictionary.
    """
    return _translate_kwargs(src, src, *keys, **aliases)


def _pop_props(src, *categories, **kwargs):
    """
    Pop out registered properties and return them in a new dictionary.
    """
    return _translate_props(src, {}, *categories, **kwargs)


def _process_props(src, *categories, **kwargs):
    """
    Translate registered properties and add translated names to the original dictionary.
    """
    return _translate_props(src, src, *categories, **kwargs)


def _pop_params(kwargs, *funcs, ignore_internal=False):
    """
    Pop parameters of the input functions or methods.
    """
    output = {}
    for func in funcs:
        sig = inspect.signature(func)
        for key in sig.parameters:
            value = kwargs.pop(key, None)
            if ignore_internal and key in INTERNAL_PARAMS:
                continue
            if value is not None:
                output[key] = value
    return output


def _fill_guide_kw(kwargs, **pairs):
    """
    Add the keyword arguments to the dictionary if not already present.
    """
    aliases = (
        ('title', 'label'),
        ('locator', 'ticks'),
        ('format', 'formatter', 'ticklabels')
    )
    for key, value in pairs.items():
        if value is None:
            continue
        keys = tuple(a for group in aliases for a in group if key in group)  # may be ()
        if not any(kwargs.get(key) is not None for key in keys):  # note any(()) is True
            kwargs[key] = value


def _guide_kw_to_arg(name, kwargs, **pairs):
    """
    Add to the `colorbar_kw` or `legend_kw` dict if there are no conflicts.
    """
    kw = kwargs.setdefault(f'{name}_kw', {})
    _fill_guide_kw(kw, **pairs)


def _guide_kw_from_obj(obj, name, kwargs):
    """
    Add to the dict from settings stored on the object if there are no conflicts.
    """
    pairs = getattr(obj, f'_{name}_kw', None)
    pairs = pairs or {}  # needed for some reason
    _fill_guide_kw(kwargs, **pairs)
    if isinstance(obj, (tuple, list, np.ndarray)):
        for iobj in obj:  # possibly iterate over matplotlib tuple/list subclasses
            _guide_kw_from_obj(iobj, name, kwargs)
    return kwargs


def _guide_kw_to_obj(obj, name, kwargs):
    """
    Add the guide keyword dict to the objects.
    """
    try:
        setattr(obj, f'_{name}_kw', kwargs)
    except AttributeError:
        pass
    if isinstance(obj, (tuple, list, np.ndarray)):
        for iobj in obj:
            _guide_kw_to_obj(iobj, name, kwargs)


class _empty_context(object):
    """
    A dummy context manager.
    """
    def __init__(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args):  # noqa: U100
        pass


class _state_context(object):
    """
    Temporarily modify attribute(s) for an arbitrary object.
    """
    def __init__(self, obj, **kwargs):
        self._obj = obj
        self._attrs_new = kwargs
        self._attrs_prev = {
            key: getattr(obj, key) for key in kwargs if hasattr(obj, key)
        }

    def __enter__(self):
        for key, value in self._attrs_new.items():
            setattr(self._obj, key, value)

    def __exit__(self, *args):  # noqa: U100
        for key in self._attrs_new.keys():
            if key in self._attrs_prev:
                setattr(self._obj, key, self._attrs_prev[key])
            else:
                delattr(self._obj, key)
