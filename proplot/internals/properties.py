#!/usr/bin/env python3
"""
Utilities for artist properties.
"""
import matplotlib.colors as mcolors
import numpy as np

from . import _not_none, docstring, warnings

# Artist property aliases. Use this rather than normalize_kwargs and _alias_maps.
# NOTE: We add the aliases 'edgewidth' and 'fillcolor' for patch edges and faces
# NOTE: Alias cannot appear as key or else _translate_kwargs will overwrite with None!
_alias_maps = {
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
    'patch': {
        'alpha': (  # secretly applies to both face and edge consistent with matplotlib
            'a', 'alphas', 'ea', 'edgealpha', 'edgealphas',
            'fa', 'facealpha', 'facealphas', 'fillalpha', 'fillalphas'
        ),
        'color': ('c', 'colors'),
        'edgecolor': ('ec', 'edgecolors'),
        'facecolor': ('fc', 'facecolors', 'fillcolor', 'fillcolors'),
        'hatch': ('h', 'hatches', 'hatching'),
        'linestyle': ('ls', 'linestyles'),
        'linewidth': ('lw', 'linewidths', 'ew', 'edgewidth', 'edgewidths'),
        'zorder': ('z', 'zorders'),
    },
    'line': {  # copied from lines.py but expanded to include plurals
        'alpha': ('a', 'alphas'),
        'color': ('c', 'colors'),
        'dashes': ('d', 'dash'),
        'drawstyle': ('ds', 'drawstyles'),
        'fillstyle': ('fs', 'fillstyles', 'mfs', 'markerfillstyle', 'markerfillstyles'),
        'linestyle': ('ls', 'linestyles'),
        'linewidth': ('lw', 'linewidths'),
        'marker': ('m', 'markers'),
        'markersize': ('s', 'ms', 'markersizes'),  # WARNING: no 'sizes' here for barb
        'markeredgewidth': ('ew', 'edgewidth', 'edgewidths', 'mew', 'markeredgewidths'),
        'markeredgecolor': ('ec', 'edgecolor', 'edgecolors', 'mec', 'markeredgecolors'),
        'markerfacecolor': (
            'fc', 'facecolor', 'facecolors', 'fillcolor', 'fillcolors',
            'mc', 'markercolor', 'markercolors', 'mfc', 'markerfacecolors'
        ),
        'zorder': ('z', 'zorders'),
    },
    'collection': {  # NOTE: collections ignore facecolor and need singular 'alpha'
        'alpha': ('a', 'alphas'),
        'colors': ('c', 'color'),
        'edgecolors': ('ec', 'edgecolor', 'mec', 'markeredgecolor', 'markeredgecolors'),
        'facecolors': (
            'fc', 'facecolor', 'fillcolor', 'fillcolors',
            'mc', 'markercolor', 'markercolors', 'mfc', 'markerfacecolor', 'markerfacecolors'  # noqa: E501
        ),
        'linestyles': ('ls', 'linestyle'),
        'linewidths': ('lw', 'linewidth', 'ew', 'edgewidth', 'edgewidths', 'mew', 'markeredgewidth', 'markeredgewidths'),  # noqa: E501
        'marker': ('m', 'markers'),
        'sizes': ('s', 'size', 'ms', 'markersize', 'markersizes'),
        'zorder': ('z', 'zorders'),
    },
    'text': {
        'color': ('c', 'fontcolor'),  # NOTE: see text.py source code
        'fontfamily': ('family', 'name', 'fontname'),
        'fontsize': ('size',),
        'fontstretch': ('stretch',),
        'fontstyle': ('style',),
        'fontvariant': ('variant',),
        'fontweight': ('weight',),
        'fontproperties': ('fp', 'font', 'font_properties'),
        'zorder': ('z', 'zorders'),
    },
}


# Unit docstrings
# NOTE: Try to fit this into a single line. Cannot break up with newline as that will
# mess up docstring indentation since this is placed in indented param lines.
_units_docstring = 'If float, units are {units}. If string, interpreted by `~proplot.utils.units`.'  # noqa: E501
docstring._snippet_manager['units.pt'] = _units_docstring.format(units='points')
docstring._snippet_manager['units.in'] = _units_docstring.format(units='inches')
docstring._snippet_manager['units.em'] = _units_docstring.format(units='em-widths')


# Artist property docstrings
# NOTE: These are needed in a few different places
_line_docstring = """
lw, linewidth, linewidths : unit-spec, default: :rc:`lines.linewidth`
    The width of the line(s).
    %(units.pt)s
ls, linestyle, linestyles : str, default: :rc:`lines.linestyle`
    The style of the line(s).
c, color, colors : color-spec, optional
    The color of the line(s). The property `cycle` is used by default.
a, alpha, alphas : float, optional
    The opacity of the line(s). Inferred from `color` by default.
"""
_patch_docstring = """
lw, linewidth, linewidths : unit-spec, default: :rc:`patch.linewidth`
    The edge width of the patch(es).
    %(units.pt)s
ls, linestyle, linestyles : str, default: '-'
    The edge style of the patch(es).
ec, edgecolor, edgecolors : color-spec, default: '{edgecolor}'
    The edge color of the patch(es).
fc, facecolor, facecolors, fillcolor, fillcolors : color-spec, optional
    The face color of the patch(es). The property `cycle` is used by default.
a, alpha, alphas : float, optional
    The opacity of the patch(es). Inferred from `facecolor` and `edgecolor` by default.
"""
_pcolor_collection_docstring = """
lw, linewidth, linewidths : unit-spec, default: 0.3
    The width of lines between grid boxes.
    %(units.pt)s
ls, linestyle, linestyles : str, default: '-'
    The style of lines between grid boxes.
ec, edgecolor, edgecolors : color-spec, default: 'k'
    The color of lines between grid boxes.
a, alpha, alphas : float, optional
    The opacity of the grid boxes. Inferred from `cmap` by default.
"""
_contour_collection_docstring = """
lw, linewidth, linewidths : unit-spec, default: 0.3 or :rc:`lines.linewidth`
    The width of the line contours. Default is ``0.3`` when adding to filled contours
    or :rc:`lines.linewidth` otherwise. %(units.pt)s
ls, linestyle, linestyles : str, default: '-' or :rc:`contour.negative_linestyle`
    The style of the line contours. Default is ``'-'`` for positive contours and
    :rcraw:`contour.negative_linestyle` for negative contours.
ec, edgecolor, edgecolors : color-spec, default: 'k' or inferred
    The color of the line contours. Default is ``'k'`` when adding to filled contours
    or inferred from `color` or `cmap` otherwise.
a, alpha, alpha : float, optional
    The opacity of the contours. Inferred from `edgecolor` by default.
"""
_text_docstring = """
name, fontname, family, fontfamily : str, optional
    The font typeface name (e.g., ``'Fira Math'``) or font family name (e.g.,
    ``'serif'``). Matplotlib falls back to the system default if not found.
size, fontsize : unit-spec or str, optional
    The font size. %(units.pt)s
    This can also be a string indicating some scaling relative to
    :rcraw:`font.size`. The sizes and scalings are shown below. The
    scalings ``'med'``, ``'med-small'``, and ``'med-large'`` are
    added by proplot while the rest are native matplotlib sizes.

    .. _font_table:

    ==========================  =====
    Size                        Scale
    ==========================  =====
    ``'xx-small'``              0.579
    ``'x-small'``               0.694
    ``'small'``, ``'smaller'``  0.833
    ``'med-small'``             0.9
    ``'med'``, ``'medium'``     1.0
    ``'med-large'``             1.1
    ``'large'``, ``'larger'``   1.2
    ``'x-large'``               1.440
    ``'xx-large'``              1.728
    ``'larger'``                1.2
    ==========================  =====

"""
docstring._snippet_manager['artist.line'] = _line_docstring
docstring._snippet_manager['artist.text'] = _text_docstring
docstring._snippet_manager['artist.patch'] = _patch_docstring.format(edgecolor='none')
docstring._snippet_manager['artist.patch_black'] = _patch_docstring.format(edgecolor='black')  # noqa: E501
docstring._snippet_manager['artist.collection_pcolor'] = _pcolor_collection_docstring
docstring._snippet_manager['artist.collection_contour'] = _contour_collection_docstring


def _get_aliases(category, *keys):
    """
    Get all available property aliases.
    """
    aliases = []
    for key in keys:
        aliases.append(key)
        aliases.extend(_alias_maps[category][key])
    return tuple(aliases)


def _list_properties(count, **kwargs):
    """
    Convert line or patch propeties to lists of values.
    """
    for key, arg in kwargs.items():
        if arg is None:
            pass
        elif (
            isinstance(arg, str)
            or not np.iterable(arg)
            or 'color' in key and mcolors.is_color_like(arg)
        ):
            arg = [arg] * count
        else:
            arg = list(arg)
            if len(arg) != count:
                raise ValueError(
                    f'Length of {key!r} properties ({len(arg)}) does '
                    f'not match the number of input arrays ({count}).'
                )
        yield arg


def _pop_properties(input, *categories, prefix=None, ignore=None, skip=None, **kwargs):
    """
    Pop the registered properties and return them in a new dictionary.
    """
    output = {}
    skip = skip or ()
    ignore = ignore or ()
    if isinstance(skip, str):  # e.g. 'sizes' for barbs() input
        skip = (skip,)
    if isinstance(ignore, str):  # e.g. 'marker' to ignore marker properties
        ignore = (ignore,)
    prefix = prefix or ''  # e.g. 'box' for boxlw, boxlinewidth, etc.
    for category in categories:
        for key, aliases in _alias_maps[category].items():
            if isinstance(aliases, str):
                aliases = (aliases,)
            opts = {
                prefix + alias: input.pop(prefix + alias, None)
                for alias in (key, *aliases)
                if alias not in skip
            }
            prop = _not_none(**opts)
            if prop is None:
                continue
            if any(string in key for string in ignore):
                warnings._warn_proplot(f'Ignoring property {key}={prop!r}.')
                continue
            if isinstance(prop, str):  # ad-hoc unit conversion
                if key in ('fontsize',):
                    from ..utils import _fontsize_to_pt
                    prop = _fontsize_to_pt(prop, **kwargs)
                if key in ('linewidth', 'linewidths', 'markersize'):
                    from ..utils import units
                    prop = units(prop, 'pt', **kwargs)
            output[key] = prop
    return output
