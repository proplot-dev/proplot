#!/usr/bin/env python3
"""
Utilities for global settings.
"""
from collections.abc import MutableMapping
from numbers import Integral, Real

import numpy as np
from cycler import Cycler
from matplotlib import RcParams
from matplotlib import rcParams as _rc_matplotlib
from matplotlib import rcParamsDefault as _rc_matplotlib_default

from . import ic  # noqa: F401
from . import _not_none, validate, warnings
from .defaults import (  # noqa: F401
    _rc_aliases,
    _rc_children,
    _rc_matplotlib_override,
    _rc_proplot_definition,
    _rc_proplot_removed,
    _rc_proplot_renamed,
)


def _get_default_param(key):
    """
    Get the default proplot rc parameter.
    """
    # NOTE: This is used for the :rc: role when compiling docs and when saving
    # proplotrc files. Includes custom proplot params and proplot overrides.
    sentinel = object()
    for dict_ in (
        _rc_proplot_default,
        _rc_matplotlib_override,  # imposed defaults
        _rc_matplotlib_default,  # native defaults
    ):
        value = dict_.get(key, sentinel)
        if value is not sentinel:
            return value
    raise KeyError(f'Invalid key {key!r}.')


def _get_param_repr(value):
    """
    Translate setting to a string suitable for saving.
    """
    # NOTE: Never safe hex strings with leading '#'. In both matplotlibrc
    # and proplotrc this will be read as comment character.
    if value is None or isinstance(value, (str, bool, Integral)):
        value = str(value)
        if value[:1] == '#':  # i.e. a HEX string
            value = value[1:]
    elif isinstance(value, Real):
        value = str(round(value, 6))  # truncate decimals
    elif isinstance(value, Cycler):
        value = repr(value)  # special case!
    elif isinstance(value, (list, tuple, np.ndarray)):
        value = ', '.join(map(_get_param_repr, value))  # sexy recursion
    else:
        value = None
    return value


def _generate_rst_table():
    """
    Return the setting names and descriptions in an RST-style table.
    """
    # Get descriptions
    descrips = {
        key: descrip for key, (_, _, descrip) in _rc_proplot_definition.items()
    }

    # Get header and divider
    keylen = 4 + len(max((*_rc_proplot_definition, 'Key'), key=len))
    vallen = len(max((*descrips.values(), 'Description'), key=len))
    spaces = 2 * ' '  # spaces between each table column
    prefix = '.. rst-class:: proplot-rctable\n\n'
    header = 'Key' + spaces + ' ' * (keylen - 3) + 'Description\n'
    divider = '=' * keylen + spaces + '=' * vallen + '\n'

    # Combine components
    string = prefix + divider + header + divider
    for key, descrip in descrips.items():
        line = '``' + key + '``' + spaces + ' ' * (keylen - len(key) - 4) + descrip
        string += line + '\n'
    string = string + divider.strip()

    return string


def _generate_yaml_section(kw, comment=True, description=False):
    """
    Return the settings as a nicely tabulated YAML-style table.
    """
    prefix = '# ' if comment else ''
    data = []
    for key, args in kw.items():
        # Possibly append description
        includes_descrip = isinstance(args, tuple) and len(args) == 3
        if not description:
            descrip = ''
            value = args[0] if includes_descrip else args
        elif includes_descrip:
            value, validator, descrip = args
            descrip = '# ' + descrip  # skip the validator
        else:
            raise ValueError(f'Unexpected input {key}={args!r}.')

        # Translate object to string
        value = _get_param_repr(value)
        if value is not None:
            data.append((key, value, descrip))
        else:
            warnings._warn_proplot(
                f'Failed to write rc setting {key} = {value!r}. Must be None, bool, '
                'string, int, float, a list or tuple thereof, or a property cycler.'
            )

    # Generate string
    string = ''
    keylen = len(max(kw, key=len))
    vallen = len(max((tup[1] for tup in data), key=len))
    for key, value, descrip in data:
        space1 = ' ' * (keylen - len(key) + 1)
        space2 = ' ' * (vallen - len(value) + 2) if descrip else ''
        string += f'{prefix}{key}:{space1}{value}{space2}{descrip}\n'
    return string.strip()


def _generate_yaml_table(changed=None, comment=None, description=False):
    """
    Return the settings as a nicely tabulated YAML-style table.
    """
    parts = [
        '#--------------------------------------------------------------------',
        '# Use this file to change the default proplot and matplotlib settings.',
        '# The syntax is identical to matplotlibrc syntax. For details see:',
        '# https://proplot.readthedocs.io/en/latest/configuration.html',
        '# https://matplotlib.org/stable/tutorials/introductory/customizing.html',
        '#--------------------------------------------------------------------',
    ]

    # User settings
    if changed is not None:  # add always-uncommented user settings
        table = _generate_yaml_section(changed, comment=False)
        parts.extend(('# Changed settings', table, ''))

    # Proplot settings
    kw = _rc_proplot_definition if description else _rc_proplot_default
    table = _generate_yaml_section(kw, description=description, comment=comment)
    parts.extend(('# Proplot settings', table, ''))

    # Matplotlib settings
    kw = _rc_matplotlib_override
    table = _generate_yaml_section(kw, comment=comment)
    parts.extend(('# Matplotlib settings', table, ''))
    return '\n'.join(parts)


def _convert_grid_param(b, key):
    """
    Translate an instruction to turn either major or minor gridlines on or off into a
    boolean and string applied to :rcraw:`axes.grid` and :rcraw:`axes.grid.which`.
    """
    ob = _rc_matplotlib['axes.grid']
    owhich = _rc_matplotlib['axes.grid.which']
    if b:
        # Gridlines are already both on, or they are off only for the
        # ones that we want to turn on. Turn on gridlines for both.
        if (
            owhich == 'both'
            or key == 'grid' and owhich == 'minor'
            or key == 'gridminor' and owhich == 'major'
        ):
            which = 'both'
        # Gridlines are off for both, or off for the ones that we
        # don't want to turn on. We can just turn on these ones.
        else:
            which = owhich
    else:
        # Gridlines are already off, or they are on for the particular
        # ones that we want to turn off. Instruct to turn both off.
        if (
            not ob
            or key == 'grid' and owhich == 'major'
            or key == 'gridminor' and owhich == 'minor'
        ):
            which = 'both'  # disable both sides
        # Gridlines are currently on for major and minor ticks, so we
        # instruct to turn on gridlines for the one we *don't* want off
        elif owhich == 'both':  # and ob is True, as already tested
            # if gridminor=False, enable major, and vice versa
            b = True
            which = 'major' if key == 'gridminor' else 'minor'
        # Gridlines are on for the ones that we *didn't* instruct to
        # turn off, and off for the ones we do want to turn off. This
        # just re-asserts the ones that are already on.
        else:
            b = True
            which = owhich
    return b, which


def _pop_settings(src, *, ignore_conflicts=True):
    """
    Pop the rc setting names and mode for a `~Configurator.context` block.
    """
    # NOTE: By default must ignore settings also present as function parameters
    # and include deprecated settings in the list.
    # NOTE: rc_mode == 2 applies only the updated params. A power user
    # could use ax.format(rc_mode=0) to re-apply all the current settings
    conflict_params = (
        'backend',
        'alpha',  # deprecated
        'color',  # deprecated
        'facecolor',  # deprecated
        'edgecolor',  # deprecated
        'linewidth',  # deprecated
        'basemap',  # deprecated
        'share',  # deprecated
        'span',  # deprecated
        'tight',  # deprecated
        'span',  # deprecated
    )
    kw = src.pop('rc_kw', None) or {}
    if 'mode' in src:
        src['rc_mode'] = src.pop('mode')
        warnings._warn_proplot(
            "Keyword 'mode' was deprecated in v0.6. Please use 'rc_mode' instead."
        )
    mode = src.pop('rc_mode', None)
    mode = _not_none(mode, 2)  # only apply updated params by default
    for key, value in tuple(src.items()):
        name = _rc_nodots.get(key, None)
        if ignore_conflicts and name in conflict_params:
            name = None  # former renamed settings
        if name is not None:
            kw[name] = src.pop(key)
    return kw, mode


class _RcParams(MutableMapping, dict):
    """
    A simple dictionary with locked inputs and validated assignments.
    """
    # NOTE: By omitting __delitem__ in MutableMapping we effectively
    # disable mutability. Also disables deleting items with pop().
    def __init__(self, source, validate):
        self._validate = validate
        for key, value in source.items():
            self.__setitem__(key, value)  # trigger validation

    def __repr__(self):
        return RcParams.__repr__(self)

    def __str__(self):
        return RcParams.__repr__(self)

    def __len__(self):
        return dict.__len__(self)

    def __iter__(self):
        # NOTE: Proplot doesn't add deprecated args to dictionary so
        # we don't have to suppress warning messages here.
        yield from sorted(dict.__iter__(self))

    def __getitem__(self, key):
        key, _ = self._check_key(key)
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        key, value = self._check_key(key, value)
        if key not in self._validate:
            raise KeyError(f'Invalid rc key {key!r}.')
        try:
            value = self._validate[key](value)
        except (ValueError, TypeError) as error:
            raise ValueError(f'Key {key!r}: {error}')
        if key is not None:
            dict.__setitem__(self, key, value)

    @staticmethod
    def _check_key(key, value=None):
        # NOTE: If we assigned from the Configurator then the deprecated key will
        # still propagate to the same 'children' as the new key.
        # NOTE: This also translates values for special cases of renamed keys.
        # Currently the special cases are 'basemap' and 'cartopy.autoextent'.
        if key in _rc_proplot_renamed:
            key_new, version = _rc_proplot_renamed[key]
            warnings._warn_proplot(
                f'The rc setting {key!r} was deprecated in version {version} and may be '  # noqa: E501
                f'removed in {warnings._next_release()}. Please use {key_new!r} instead.'  # noqa: E501
            )
            if key == 'basemap':  # special case
                value = ('cartopy', 'basemap')[int(bool(value))]
            if key == 'cartopy.autoextent':
                value = ('globe', 'auto')[int(bool(value))]
            key = key_new
        if key in _rc_proplot_removed:
            info, version = _rc_proplot_removed[key]
            raise KeyError(
                f'The rc setting {key!r} was removed in version {version}.'
                + (info and ' ' + info)
            )
        return key, value

    def copy(self):
        source = {key: dict.__getitem__(self, key) for key in self}
        return _RcParams(source, self._validate)


# Validate the default settings dictionaries using a custom proplot _RcParams and the
# original matplotlib RcParams. Also surreptitiously add proplot font settings to the
# font keys list (beoolean below always evalutes to True) font keys list during init.
_rc_proplot_default = {
    key: value
    for key, (value, _, _) in _rc_proplot_definition.items()
}
_rc_proplot_validate = {
    key: validator
    for key, (_, validator, _) in _rc_proplot_definition.items()
    if not (validator is validate._validate_fontsize and validate.FONT_KEYS.add(key))
}
_rc_proplot_default = _RcParams(_rc_proplot_default, _rc_proplot_validate)
_rc_matplotlib_override = RcParams(_rc_matplotlib_override)

# Important joint matplotlib proplot constants
# NOTE: The 'nodots' dictionary should include removed and renamed settings
_rc_categories = {
    '.'.join(name.split('.')[:i + 1])
    for dict_ in (
        _rc_proplot_default,
        _rc_matplotlib_default
    )
    for name in dict_
    for i in range(len(name.split('.')) - 1)
}
_rc_nodots = {
    name.replace('.', ''): name
    for dict_ in (
        _rc_proplot_definition,
        _rc_proplot_renamed,
        _rc_proplot_removed,
        _rc_matplotlib_default,
    )
    for name in dict_.keys()
}
