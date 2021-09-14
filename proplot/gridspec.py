#!/usr/bin/env python3
"""
The gridspec and subplot grid classes used throughout ProPlot.
"""
import functools
import inspect
import itertools
import re
from collections.abc import MutableSequence
from numbers import Integral

import matplotlib.axes as maxes
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import numpy as np

from . import axes as paxes
from .config import rc
from .internals import ic  # noqa: F401
from .internals import _not_none, docstring, warnings
from .utils import _fontsize_to_pt, units

__all__ = [
    'GridSpec',
    'SubplotGrid',
    'SubplotsContainer'  # deprecated
]


# Gridspec vector arguments
# Valid for figure() and GridSpec()
_shared_docstring = """
left, right, top, bottom : unit-spec, optional
    The fixed space between the subplots and the figure edge. Default is ``None``.
    %(units.em)s
    If ``None``, the space is determined automatically based on the font size and
    axis sharing settings. If :rcraw:`subplots.tight` is ``True``, the space is
    determined by the tight layout algorithm.
"""
_scalar_docstring = """
wspace, hspace, space : unit-spec, optional
    The fixed space between grid columns, rows, or both. Default is ``None``.
    %(units.em)s
    If ``None``, the space is determined automatically based on the font size and
    axis sharing settings. If :rcraw:`subplots.tight` is ``True``, the space is
    determined by the tight layout algorithm.
"""
_vector_docstring = """
wspace, hspace, space : unit-spec or sequence, optional
    The fixed space between grid columns, rows, and both, respectively. If
    float, string, or ``None``, this value is expanded into lists of length
    ``ncols - 1`` (for `wspace`) or length ``nrows - 1`` (for `hspace`). If
    a sequence, its length must match these lengths. Default is ``None``.
    %(units.em)s

    For elements equal to ``None``, the space is determined automatically based
    on the font and tick settings. If :rcraw:`subplots.tight` is ``True``, the
    space is determined by the tight layout algorithm. Otherwise, a sensible default
    value is chosen. For example, ``subplots(ncols=3, tight=True, wspace=(2, None))``
    fixes the space between columns 1 and 2 but lets the tight layout algorithm
    determine the space between columns 2 and 3.
wratios, hratios : float or sequence, optional
    Passed to `~proplot.gridspec.GridSpec`, denotes the width and height
    ratios for the subplot grid. Length of `wratios` must match the number
    of rows, and length of `hratios` must match the number of columns.
width_ratios, height_ratios
    Aliases for `wratios`, `hratios`. Included for consistency with
    the `matplotlib.pyplot.subplots` command.
wpad, hpad, pad : unit-spec or sequence, optional
    The tight layout padding between columns, rows, and both, respectively. Unlike
    ``space``, these control the padding between subplot content (including text,
    ticks, etc.) rather than subplot edges. As with ``space``, these can be scalars
    or arrays optionally containing ``None``. Default is `innerpad`.
    %(units.em)s
"""
_tight_docstring = """
wequal, hequal, equal :  bool, optional
    Whether to make the tight layout algorithm apply equal spacing between columns,
    rows, or both. Default is ``False``. Ignored if :rcraw:`tight` is ``False``.
outerpad : unit-spec, optional
    The tight layout padding around the left, right, top, and bottom edges
    of the figure.  Default is :rc:`subplots.outerpad`.
    %(units.em)s
innerpad : unit-spec, optional
    The scalar tight layout padding between columns and rows. Synonymous with
    `pad`. Default is :rc:`subplots.innerpad`.
    %(units.em)s
panelpad : unit-spec, optional
    The tight layout padding between subplots and axes panels and between "stacked"
    panels. Default is :rc:`subplots.panelpad`.
    %(units.em)s
"""
docstring._snippet_manager['gridspec.shared'] = _shared_docstring
docstring._snippet_manager['gridspec.scalar'] = _scalar_docstring
docstring._snippet_manager['gridspec.vector'] = _vector_docstring
docstring._snippet_manager['gridspec.tight'] = _tight_docstring


def _disable_method(attr):
    """
    Disable the inherited method.
    """
    def _dummy_method(*args):
        raise AttributeError(f'Method {attr}() is disabled on ProPlot gridspecs.')
    _dummy_method.__name__ = attr
    return _dummy_method


class _SubplotSpec(mgridspec.SubplotSpec):
    """
    A thin `~matplotlib.gridspec.SubplotSpec` subclass with a nice string
    representation and a few helper methods.
    """
    def __repr__(self):
        # NOTE: Also include panel obfuscation here to avoid confusion. If this
        # is a panel slot generated internally then show zero info.
        try:
            nrows, ncols, num1, num2 = self._get_subplot_geometry()
        except (IndexError, ValueError, AttributeError):
            return 'SubplotSpec(unknown)'
        else:
            return f'SubplotSpec(nrows={nrows}, ncols={ncols}, index=({num1}, {num2}))'

    def _get_geometry(self):
        """
        Return the full geometry.
        """
        return super().get_geometry()

    def _get_subplot_geometry(self):
        """
        Return the subplot geometry. May trigger an error if this is a panel slot.
        """
        gs = self.get_gridspec()
        nrows, ncols = gs.get_subplot_geometry()
        num1, num2 = gs._convert_full_to_subplot(self.num1, self.num2)
        return nrows, ncols, num1, num2

    def _get_rows_columns(self, ncols=None):
        """
        Return the full rows and columns.
        """
        _, ncols_true, num1, num2 = self._get_geometry()
        ncols = _not_none(ncols, ncols_true)
        row1, col1 = divmod(num1, ncols)
        row2, col2 = divmod(num2, ncols)
        return row1, row2, col1, col2

    def _get_subplot_rows_columns(self, ncols=None):
        """
        Return the subplot rows and columns. May trigger error if this is a panel slot.
        """
        _, ncols_true, num1, num2 = self._get_subplot_geometry()
        ncols = _not_none(ncols, ncols_true)
        row1, col1 = divmod(num1, ncols)
        row2, col2 = divmod(num2, ncols)
        return row1, row2, col1, col2

    def get_position(self, *args, **kwargs):
        # Silent override. Older matplotlib versions can create subplots
        # with negative heights and widths that crash on instantiation.
        # Instead better to dynamically adjust the bounding box and hope
        # that subsequent adjustments will correct the subplot position.
        result = super().get_position(*args, **kwargs)
        bbox = result[0] if isinstance(result, tuple) else result
        if isinstance(bbox, mtransforms.BboxBase):
            extents = bbox.extents
            extents = [
                extents[0],
                extents[1],
                max(extents[2], extents[0]),
                max(extents[3], extents[1]),
            ]
            bbox = mtransforms.Bbox.from_extents(extents)
            if isinstance(result, tuple):
                result = (bbox, *result[1:])
            else:
                result = bbox
        return result


class GridSpec(mgridspec.GridSpec):
    """
    A `~matplotlib.gridspec.GridSpec` subclass that permits variable spacing
    between successive rows and columns and hides "panel slots" from indexing.
    """
    def __repr__(self):
        srows, scols = self.get_subplot_geometry()
        prows, pcols = self.get_panel_geometry()
        params = {'nrows': srows, 'ncols': scols}
        if prows:
            params['panelrows'] = prows
        if pcols:
            params['panelcols'] = pcols
        params = ', '.join(f'{key}={value!r}' for key, value in params.items())
        return f'GridSpec({params})'

    def __getattr__(self, attr):
        # Redirect to private 'layout' attributes that are fragile w.r.t.
        # matplotlib version. Cannot set these by calling super().__init__()
        # because we make spacing arguments non-settable properties.
        if 'layout' in attr:
            return None
        super().__getattribute__(attr)  # native error message

    @docstring._snippet_manager
    def __init__(self, nrows=1, ncols=1, **kwargs):
        """
        Parameters
        ----------
        nrows : int, optional
            The number of rows in the subplot grid.
        ncols : int, optional
            The number of columns in the subplot grid.

        Other parameters
        ----------------
        %(gridspec.shared)s
        %(gridspec.vector)s
        %(gridspec.tight)s

        Note
        ----
        Adding axes panels, axes or figure colorbars, and axes or figure legends
        quietly augments the gridspec geometry by inserting "panel slots". However
        subsequently indexing the gridspec with ``gs[num]`` or ``gs[row, col]`` will
        ignore the "panel slots". This permits adding new subplots by passing
        ``gs[num]`` or ``gs[row, col]`` to `~proplot.figure.Figure.add_subplot`
        even in the presence of panels. See `~GridSpec.__getitem__` for details.
        """
        # Gridspec properties
        self._nrows = nrows
        self._ncols = ncols
        self._left = None
        self._right = None
        self._bottom = None
        self._top = None
        self._hspace = [None] * (nrows - 1)
        self._wspace = [None] * (ncols - 1)
        self._hratios = [1] * nrows
        self._wratios = [1] * ncols
        self._left_default = None
        self._right_default = None
        self._bottom_default = None
        self._top_default = None
        self._hspace_default = [None] * (nrows - 1)
        self._wspace_default = [None] * (ncols - 1)
        self._figure = None  # initial state

        # Capture rc settings used for default spacing
        # NOTE: This is consistent with conversion of 'em' units to inches on gridspec
        # instantiation. In general it seems strange for future changes to rc settings
        # to magically update an existing gridspec layout. This also may improve
        # lookup time since get_grid_positions() is called heaviliy.
        scales = {'in': 0, 'inout': 0.5, 'out': 1, None: 1}
        self._xtickspace = scales[rc['xtick.direction']] * rc['xtick.major.size']
        self._ytickspace = scales[rc['ytick.direction']] * rc['ytick.major.size']
        self._xticklabelspace = _fontsize_to_pt(rc['xtick.labelsize']) + rc['xtick.major.pad']  # noqa: E501
        self._yticklabelspace = 3 * _fontsize_to_pt(rc['ytick.labelsize']) + rc['ytick.major.pad']  # noqa: E501
        self._labelspace = _fontsize_to_pt(rc['axes.labelsize']) + rc['axes.labelpad']
        self._titlespace = _fontsize_to_pt(rc['axes.titlesize']) + rc['axes.titlepad']

        # Tight layout and panel-related properties
        # NOTE: The wpanels and hpanels contain empty strings '' (indicating main axes),
        # or one of 'l', 'r', 'b', 't' (indicating axes panels) or 'f' (figure panels)
        outerpad = _not_none(kwargs.pop('outerpad', None), rc['subplots.outerpad'])
        innerpad = _not_none(kwargs.pop('innerpad', None), rc['subplots.innerpad'])
        panelpad = _not_none(kwargs.pop('panelpad', None), rc['subplots.panelpad'])
        pad = _not_none(kwargs.pop('pad', None), innerpad)  # alias of innerpad
        self._outerpad = units(outerpad, 'em', 'in')
        self._innerpad = units(innerpad, 'em', 'in')
        self._panelpad = units(panelpad, 'em', 'in')
        self._hpad = [units(pad, 'em', 'in')] * (nrows - 1)
        self._wpad = [units(pad, 'em', 'in')] * (ncols - 1)
        self._hequal = False
        self._wequal = False
        self._hpanels = [''] * nrows  # axes and figure panel identification
        self._wpanels = [''] * ncols
        self._fpanels = {  # array representation of figure panel spans
            'left': np.empty((0, nrows), dtype=bool),
            'right': np.empty((0, nrows), dtype=bool),
            'bottom': np.empty((0, ncols), dtype=bool),
            'top': np.empty((0, ncols), dtype=bool),
        }
        self._update_params(pad=pad, **kwargs)

    def __getitem__(self, key):
        """
        Get a `~matplotlib.gridspec.SubplotSpec`. Slots allocated for axes panels,
        colorbars, and legends are ignored. For example, given a gridspec with 3
        subplot rows, 3 subplot columns, 6 panel rows, and 6 panel columns, calling
        ``gs[1, 1]`` returns a `~matplotlib.gridspec.SubplotSpec` corresponding to
        the second *subplot* row and column rather than a *panel* row or column.
        """
        return self._get_subplot_spec(key, includepanels=False)

    def _convert_full_to_subplot(self, *args, which=None):
        """
        Convert indices from the full geometry into the "subplot" gridspec
        geometry. If `which` is not passed these should be flattened indices.
        """
        nums = []
        idxs = self._get_subplot_indices(which)
        for arg in args:
            try:
                nums.append(idxs.index(arg))
            except ValueError:
                raise ValueError(f'Invalid gridspec index {arg}.')
        return nums[0] if len(nums) == 1 else nums

    def _convert_subplot_to_full(self, *args, which=None):
        """
        Convert indices from the "subplot" gridspec geometry into indices for the full
        geometry. If `which` is not passed these should be flattened indices.
        """
        nums = []
        idxs = self._get_subplot_indices(which)
        for arg in args:
            try:
                nums.append(idxs[arg])
            except (IndexError, TypeError):
                raise ValueError(f'Invalid gridspec index {arg}.')
        return nums[0] if len(nums) == 1 else nums

    def _get_current_space(self, key):
        """
        Get the currently active space accounting for both default
        values and explicit user-specified values.
        """
        # NOTE: Default panel spaces should already have been filled by _insert_panel.
        # They use 'panelpad' and the panel-local 'share' setting. This function
        # instead fills spaces between subplots depending on sharing setting.
        fig = self.figure
        if not fig:
            raise ValueError('Figure must be assigned to get grid positions.')
        attr = f'_{key}'  # user-specified
        attr_default = f'_{key}_default'  # default values
        value = getattr(self, attr)
        value_default = getattr(self, attr_default)
        if key in ('left', 'right', 'bottom', 'top'):
            if value_default is None:
                value_default = self._get_default_space(key)
                setattr(self, attr_default, value_default)
            return _not_none(value, value_default)
        elif key in ('wspace', 'hspace'):
            result = []
            for i, (val, val_default) in enumerate(zip(value, value_default)):
                if val_default is None:
                    val_default = self._get_default_space(key)
                    value_default[i] = val_default
                result.append(_not_none(val, val_default))
            return result
        else:
            raise ValueError(f'Unknown space parameter {key!r}.')

    def _get_default_space(self, key, pad=None, share=None, title=True):
        """
        Return suitable default spacing given a shared axes setting.
        This is only relevant when "tight layout" is disabled.
        """
        # NOTE: Internal spacing args are stored in inches to simplify the
        # get_grid_positions() calculations.
        fig = self.figure
        if fig is None:
            raise RuntimeError('Figure must be assigned.')
        if key == 'right':
            pad = _not_none(pad, self._outerpad)
            space = 0
        elif key == 'top':
            pad = _not_none(pad, self._outerpad)
            space = self._titlespace if title else 0
        elif key == 'left':
            pad = _not_none(pad, self._outerpad)
            space = self._labelspace + self._yticklabelspace + self._ytickspace
        elif key == 'bottom':
            pad = _not_none(pad, self._outerpad)
            space = self._labelspace + self._xticklabelspace + self._xtickspace
        elif key == 'wspace':
            pad = _not_none(pad, self._innerpad)
            share = _not_none(share, fig._sharey, 0)
            space = self._ytickspace
            if share < 3:
                space += self._yticklabelspace
            if share < 1:
                space += self._labelspace
        elif key == 'hspace':
            pad = _not_none(pad, self._innerpad)
            share = _not_none(share, fig._sharex, 0)
            space = self._xtickspace
            if title:
                space += self._titlespace
            if share < 3:
                space += self._xticklabelspace
            if share < 1:
                space += self._labelspace
        else:
            raise ValueError(f'Invalid space key {key!r}.')
        return pad + space / 72

    def _get_panel_list(self, which=None):
        """
        Get the panel list for the rows, columns, or flattened array.
        """
        hpanels = self._hpanels
        wpanels = self._wpanels
        if which is None:
            panels = tuple(h + w for h, w in itertools.product(hpanels, wpanels))
        elif which in 'xw':
            panels = wpanels
        elif which in 'yh':
            panels = hpanels
        else:
            raise ValueError(f'Invalid which={which!r}.')
        return panels

    def _get_panel_indices(self, which=None, space=False):
        """
        Get the indices associated with "panel" gridspec slots or spaces.
        """
        panels = self._get_panel_list(which)
        if not space:
            idxs = [i for i, s in enumerate(panels) if s]
        else:
            idxs = [
                i for i, (p1, p2) in enumerate(zip(panels[:-1], panels[1:]))
                if p1 == 'f' and p2 == 'f'
                or p1 in ('l', 't') and p2 in ('l', 't', '')
                or p1 in ('r', 'b', '') and p2 in ('r', 'b')
            ]
        return idxs

    def _get_subplot_indices(self, which=None, space=False):
        """
        Get the indices associated with "subplot" gridspec slots or spaces.
        """
        panels = self._get_panel_list(which)
        length = len(panels) - 1 if space else len(panels)
        idxs = self._get_panel_indices(which=which, space=space)
        idxs = [i for i in range(length) if i not in idxs]
        return idxs

    def _get_subplot_spec(self, key, includepanels=False):
        """
        Generate a subplotspec either ignoring panels or including panels.
        """
        # Convert the indices into endpoint-inclusive (start, stop)
        def _normalize_index(key, size, axis=None):  # noqa: E306
            if isinstance(key, slice):
                start, stop, _ = key.indices(size)
                if stop > start:
                    return start, stop - 1
            else:
                if key < 0:
                    key += size
                if 0 <= key < size:
                    return key, key  # endpoing inclusive
            extra = 'for gridspec' if axis is None else f'along axis {axis}'
            raise IndexError(f'Invalid index {key} {extra} with size {size}.')

        # Normalize the indices
        if includepanels:
            nrows, ncols = self.get_geometry()
        else:
            nrows, ncols = self.get_subplot_geometry()
        if not isinstance(key, tuple):  # usage gridspec[1,2]
            num1, num2 = _normalize_index(key, nrows * ncols)
        elif len(key) == 2:
            k1, k2 = key
            num1 = _normalize_index(k1, nrows, axis=0)
            num2 = _normalize_index(k2, ncols, axis=1)
            num1, num2 = np.ravel_multi_index((num1, num2), (nrows, ncols))
        else:
            raise ValueError(f'Invalid index {key!r}.')

        # Return the subplotspec
        if not includepanels:
            num1, num2 = self._convert_subplot_to_full(num1, num2)
        return _SubplotSpec(self, num1, num2)

    def _parse_axes_panel(self, side, ax):
        """
        Return the indices associated with a new axes panel on the specified side.
        """
        # Get gridspec and subplotspec indices
        ss = ax.get_subplotspec()
        offset = len(ax._panel_dict[side]) + 1
        row1, row2, col1, col2 = ss._get_rows_columns()
        if side in ('left', 'right'):
            iratio = col1 - offset if side == 'left' else col2 + offset
            start, stop = row1, row2 + 1
        else:
            iratio = row1 - offset if side == 'top' else row2 + offset
            start, stop = col1, col2 + 1

        # Return subplotspec indices
        return iratio, slice(start, stop)

    def _parse_figure_panel(self, side, span):
        """
        Return the indices associated with a new figure panel on the specified side.
        Try to find room in the current mosaic of figure panels.
        """
        # Parse the input span
        # NOTE: Here the '_fpanels' array always retains the original gridspec
        # size before panels were drawn.
        # NOTE: Here the 'span' indices start at '1' by analogy with add_subplot()
        # integers and with main subplot numbers. Also *ignores panel slots*.
        array = self._fpanels[side]
        nacross = self._ncols if side in ('left', 'right') else self._nrows
        npanels, nalong = array.shape
        span = _not_none(span, (1, nalong))
        span = np.atleast_1d(span)
        if span.size not in (1, 2):
            raise ValueError(f'Invalid span={span!r}. Must be scalar or 2-tuple of coordinates.')  # noqa: E501
        if any(s < 1 or s > nalong for s in span):
            raise ValueError(f'Invalid span={span!r}. Coordinates must satisfy 1 <= c <= {nalong}.')  # noqa: E501

        # Modify the array
        # NOTE: The array is array of booleans, where each row corresponds to a figure
        # panel, moving toward the outside, and True indicates the slot is filled.
        start, stop = span[0] - 1, span[-1]  # non-inclusive starting at zero
        iratio = -1 if side in ('left', 'top') else nacross  # default values
        for i in range(npanels):
            if any(array[i, start:stop]):  # filled
                continue
            array[i, start:stop] = True
            if side in ('left', 'top'):  # descending moves us closer to 0
                iratio = npanels - 1 - i  # index in ratios array
            else:  # descending array moves us closer to nacross - 1
                iratio = nacross - (npanels - i)  # index in ratios array
            break
        if iratio == -1 or iratio == nacross:  # no slots found so we must add to array
            iarray = np.zeros((1, nalong), dtype=bool)
            iarray[0, start:stop] = True
            array = np.concatenate((array, iarray), axis=0)
            self._fpanels[side] = array  # update array

        # Return subplotspec indices
        # NOTE: Convert using the lengthwise indices
        which = 'h' if side in ('left', 'right') else 'w'
        start, stop = self._convert_subplot_to_full(start, stop - 1, which=which)
        return iratio, slice(start, stop + 1)

    def _insert_panel(
        self, side, arg, *,
        share=None, width=None, space=None, pad=None, filled=False,
    ):
        """
        Insert a panel slot into the existing gridspec. The `side` is the panel side
        and the `arg` is either an axes instance or the figure row-column span.
        Subsequently indexing the gridspec will ignore the slots occupied by panels.
        """
        fig = self.figure
        if fig is None:
            raise RuntimeError('Figure must be assigned to gridspec.')
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid side {side}.')
        if side in ('left', 'right'):
            pads = self._wpad
            panels = self._wpanels
            ratios = self._wratios
            spaces = self._wspace
            spaces_default = self._wspace_default
        else:
            pads = self._hpad
            panels = self._hpanels
            ratios = self._hratios
            spaces = self._hspace
            spaces_default = self._hspace_default

        # Get the subplotspec index for the panel
        if isinstance(arg, maxes.SubplotBase):
            slot_type = side[0]
            idx, span = self._parse_axes_panel(side, arg)
        else:
            slot_type = 'f'
            idx, span = self._parse_figure_panel(side, arg)
        idx_off = 1 * bool(side in ('top', 'left'))
        idx_space = idx - 1 * bool(side in ('bottom', 'right'))

        # Get user-input properties and default properties
        # NOTE: For panels there are no sharing 'levels' just boolean toggle. Also
        # note this gives totally wrong space for 'right' and 'top' colorbars with
        # tick locations on the outside but that's ok... otherwise would have to
        # modify existing space which is overkill. People should use tight layout
        pad = units(pad, 'em', 'in')
        space = units(space, 'em', 'in')
        width = units(width, 'in')
        share = False if filled else _not_none(share, True)
        pad_default = (
            self._panelpad
            if slot_type != 'f'
            or side in ('left', 'top') and panels[0] == 'f'
            or side in ('right', 'bottom') and panels[-1] == 'f'
            else self._innerpad
        )
        space_default = (
            _not_none(pad, pad_default)
            if side in ('top', 'right')
            else self._get_default_space(
                'hspace' if side == 'bottom' else 'wspace',
                title=False,  # no title space
                share=3 if share else 0,
                pad=_not_none(pad, pad_default),
            )
        )
        width_default = units(
            rc['colorbar.width' if filled else 'subplots.panelwidth'], 'in'
        )

        # Adjust space, ratio, and panel indicator arrays
        ncols = self.ncols
        newrow = newcol = None
        slot_exists = idx not in (-1, len(panels)) and panels[idx] == slot_type
        if slot_exists:
            # Slot already exists
            # Overwrite width, pad, or space only if they were provided by the user
            spaces_default[idx_space] = space_default
            if width is not None:
                ratios[idx] = width
            if pad is not None:
                pads[idx_space] = pad
            if space is not None:
                spaces[idx_space] = space
        else:
            # Modify basic geometry
            idx += idx_off
            idx_space += idx_off
            if side in ('left', 'right'):
                newcol = idx
                self._ncols += 1
            else:
                newrow = idx
                self._nrows += 1
            panels.insert(idx, slot_type)
            ratios.insert(idx, _not_none(width, width_default))
            pads.insert(idx_space, _not_none(pad, pad_default))
            spaces.insert(idx_space, space)
            spaces_default.insert(idx_space, space_default)

        # Update the SubplotSpec occuped by the axes
        inserts = (newrow, newrow, newcol, newcol)
        for ax in fig._iter_axes(hidden=True, children=True):
            # Get old index
            # NOTE: Endpoints are inclusive, not exclusive!
            if not isinstance(ax, maxes.SubplotBase):
                continue
            gs = ax.get_subplotspec().get_gridspec()
            ss = ax.get_subplotspec().get_topmost_subplotspec()
            # Get a new subplotspec
            coords = list(ss._get_rows_columns(ncols))
            for i in range(4):
                if inserts[i] is not None and coords[i] >= inserts[i]:
                    coords[i] += 1
            row1, row2, col1, col2 = coords
            key1 = slice(row1, row2 + 1)
            key2 = slice(col1, col2 + 1)
            ss_new = self._get_subplot_spec((key1, key2), includepanels=True)
            # Apply new subplotspec
            # NOTE: We should only have one possible level of GridSpecFromSubplotSpec
            # nesting -- from making side colorbars with length less than 1.
            if ss is ax.get_subplotspec():
                ax.set_subplotspec(ss_new)
            elif ss is getattr(gs, '_subplot_spec', None):
                gs._subplot_spec = ss_new
            else:
                raise RuntimeError('Unexpected GridSpecFromSubplotSpec nesting.')
            ax._reposition_subplot()

        # Update the figure size and layout
        figsize = self._calc_figsize()
        if figsize is not None:
            fig.set_size_inches(figsize, internal=True, forward=False)
        else:
            self.update()

        # Return a subplotspec
        # NOTE: For figure panels indices are determined by user-input spans.
        key = (span, idx) if side in ('left', 'right') else (idx, span)
        ss = self._get_subplot_spec(key, includepanels=True)  # bypass panel obfuscation
        return ss, share

    def _calc_figsize(self):
        """
        Return an updated auto layout figure size accounting for
        gridspec and figure parameters
        """
        fig = self.figure
        if fig is None:  # drawing before subplots are added?
            return
        ax = fig._subplot_dict.get(fig._refnum, None)
        if ax is None:  # drawing before subplots are added?
            return
        y1, y2, x1, x2 = ax.get_subplotspec()._get_rows_columns()
        refhspace = sum(self.hspace[y1:y2])
        refwspace = sum(self.wspace[x1:x2])
        refhpanel = sum(self.hratios[i] for i in range(y1, y2 + 1) if self._hpanels[i])
        refwpanel = sum(self.wratios[i] for i in range(x1, x2 + 1) if self._wpanels[i])
        refhsubplot = sum(self.hratios[i] for i in range(y1, y2 + 1) if not self._hpanels[i])  # noqa: E501
        refwsubplot = sum(self.wratios[i] for i in range(x1, x2 + 1) if not self._wpanels[i])  # noqa: E501

        # Get the reference sizes
        # NOTE: The sizing arguments should have been normalized already
        figwidth, figheight = fig._figwidth, fig._figheight
        refwidth, refheight = fig._refwidth, fig._refheight
        refaspect = _not_none(fig._refaspect, fig._refaspect_default)
        if refheight is None and figheight is None:
            if figwidth is not None:
                gridwidth = figwidth - self.spacewidth - self.panelwidth
                refwidth = gridwidth * refwsubplot / self.subplotwidth
            if refwidth is not None:  # WARNING: do not change to elif!
                refheight = refwidth / refaspect
            else:
                raise RuntimeError('Figure size arguments are all missing.')
        if refwidth is None and figwidth is None:
            if figheight is not None:
                gridheight = figheight - self.spaceheight - self.panelheight
                refheight = gridheight * refhsubplot / self.subplotheight
            if refheight is not None:
                refwidth = refheight * refaspect
            else:
                raise RuntimeError('Figure size arguments are all missing.')

        # Get the auto figure size. Might trigger 'not enough room' error later
        # NOTE: For e.g. [[1, 1, 2, 2], [0, 3, 3, 0]] we make sure to still scale the
        # reference axes like a square even though takes two columns of gridspec.
        if refheight is not None:
            refheight -= refhspace + refhpanel
            gridheight = refheight * self.subplotheight / refhsubplot
            figheight = gridheight + self.spaceheight + self.panelheight
        if refwidth is not None:
            refwidth -= refwspace + refwpanel
            gridwidth = refwidth * self.subplotwidth / refwsubplot
            figwidth = gridwidth + self.spacewidth + self.panelwidth

        # Return the figure size
        figsize = (figwidth, figheight)
        if all(np.isfinite(figsize)):
            return figsize
        else:
            warnings._warn_proplot(f'Auto resize failed. Invalid figsize {figsize}.')

    def _calc_space(self, w):
        """
        Get tight layout spaces between the input subplot rows or columns.
        """
        # Get constants
        fig = self.figure
        if not fig:
            return
        if w == 'w':
            x, y = 'xy'
            nacross = self.nrows
            space = self.wspace
            pad = self.wpad
        else:
            x, y = 'yx'
            nacross = self.ncols
            space = self.hspace
            pad = self.hpad

        # Iterate along each row or column space
        axs = tuple(fig._iter_axes(hidden=True, children=False))
        space = list(space)  # a copy
        ralong = np.array([ax._range_subplotspec(x) for ax in axs])
        racross = np.array([ax._range_subplotspec(y) for ax in axs])
        for i, (s, p) in enumerate(zip(space, pad)):
            # Find axes that abutt aginst this row or column space
            groups = []
            filt1 = ralong[:, 1] == i  # i.e. r / b edge abutts against this
            filt2 = ralong[:, 0] == i + 1  # i.e. l / t edge abutts against this
            for j in range(nacross):  # e.g. each row
                # Get the indices for axes that meet this row or column edge.
                filt = (racross[:, 0] <= j) & (j <= racross[:, 1])
                if sum(filt) < 2:
                    continue  # no interface
                idx1, = np.where(filt & filt1)
                idx2, = np.where(filt & filt2)
                if idx1.size != 1 or idx2.size != 1:
                    continue
                idx1, idx2 = idx1[0], idx2[0]
                # Put axes into unique groups and store as (l, r) or (b, t) pairs.
                ax1, ax2 = axs[idx1], axs[idx2]
                if x != 'x':  # order bottom-to-top
                    ax1, ax2 = ax2, ax1
                newgroup = True
                for (group1, group2) in groups:
                    if ax1 in group1 or ax2 in group2:
                        newgroup = False
                        group1.add(ax1)
                        group2.add(ax2)
                        break
                if newgroup:
                    groups.append([{ax1}, {ax2}])  # form new group
            # Determing the spaces using cached tight bounding boxes
            # NOTE: Set gridspec space to zero if there are no adjacent edges
            margins = []
            for (group1, group2) in groups:
                x1 = max(ax._range_tightbbox(x)[1] for ax in group1)
                x2 = min(ax._range_tightbbox(x)[0] for ax in group2)
                margins.append((x2 - x1) / self.figure.dpi)
            s = 0 if not margins else max(0, s - min(margins) + p)
            space[i] = s

        return space

    def _auto_layout_aspect(self):
        """
        Update the underlying default aspect ratio.
        """
        # Get the axes
        fig = self.figure
        if not fig:
            return
        ax = fig._subplot_dict.get(fig._refnum, None)
        if ax is None:
            return

        # Get aspect ratio
        ratio = ax.get_aspect()  # the aspect ratio in *data units*
        if ratio == 'auto':
            return
        elif ratio == 'equal':
            ratio = 1
        elif isinstance(ratio, str):
            raise RuntimeError(f'Unknown aspect ratio mode {ratio!r}.')

        # Compare to current aspect
        xscale, yscale = ax.get_xscale(), ax.get_yscale()
        if xscale == 'linear' and yscale == 'linear':
            aspect = ratio / ax.get_data_ratio()
        elif xscale == 'log' and yscale == 'log':
            aspect = ratio / ax.get_data_ratio_log()
        else:
            return  # matplotlib should have issued warning
        if fig._refaspect is not None:
            return  # fixed by user
        if np.isclose(aspect, fig._refaspect_default):
            return  # close enough to the default aspect
        fig._refaspect_default = aspect

        # Update the layout
        figsize = self._calc_figsize()
        if not fig._is_same_size(figsize):
            fig.set_size_inches(figsize, internal=True)

    def _auto_layout_space(self, renderer):
        """
        Update the underlying spaces with tight layout values. If `resize` is
        ``True`` and the auto figure size has changed then update the figure
        size. Either way always update the subplot positions.
        """
        # Initial stuff
        fig = self.figure
        if not fig:
            return
        if not any(fig._iter_axes(hidden=True, children=False)):
            return  # skip tight layout if there are no subplots in the figure

        # Get the tight bounding box around the whole figure.
        # NOTE: This triggers proplot.axes.Axes.get_tightbbox which *caches* the
        # computed bounding boxes used by _range_tightbbox below.
        pad = self._outerpad
        obox = fig.bbox_inches  # original bbox
        bbox = fig.get_tightbbox(renderer)

        # Calculate new figure margins
        # NOTE: Negative spaces are common where entire rows/columns of gridspec
        # are empty but it seems to result in wrong figure size + grid positions. Not
        # worth correcting so instead enforce positive margin sizes. Will leave big
        # empty slot but that is probably what should happen under this scenario.
        left = self.left
        bottom = self.bottom
        right = self.right
        top = self.top
        self._left_default = max(0, left - (bbox.xmin - 0) + pad)
        self._bottom_default = max(0, bottom - (bbox.ymin - 0) + pad)
        self._right_default = max(0, right - (obox.xmax - bbox.xmax) + pad)
        self._top_default = max(0, top - (obox.ymax - bbox.ymax) + pad)

        # Calculate new subplot row and column spaces. Enforce equal
        # default spaces between main subplot edges if requested.
        hspace = self._calc_space('h')
        wspace = self._calc_space('w')
        if self._hequal:
            idxs = self._get_subplot_indices('h', space=True)
            space = max(hspace[i] for i in idxs)
            for i in idxs:
                hspace[i] = space
        if self._wequal:
            idxs = self._get_subplot_indices('w', space=True)
            space = max(wspace[i] for i in idxs)
            for i in idxs:
                wspace[i] = space
        self._hspace_default = hspace
        self._wspace_default = wspace

        # Update the layout
        # NOTE: fig.set_size_inches() always updates the gridspec to enforce fixed
        # spaces (necessary since native position coordinates are figure-relative)
        # and to enforce fixed panel ratios. So only self.update() if we skip resize.
        figsize = self._calc_figsize()
        if not fig._is_same_size(figsize):
            fig.set_size_inches(figsize, internal=True)
        else:
            self.update()

    def _update_params(
        self, *,
        left=None, bottom=None, right=None, top=None,
        wspace=None, hspace=None, space=None,
        wpad=None, hpad=None, pad=None,
        wequal=None, hequal=None, equal=None,
        outerpad=None, innerpad=None, panelpad=None,
        hratios=None, wratios=None, width_ratios=None, height_ratios=None,
    ):
        """
        Update the user-specified properties.
        """
        # Assign scalar args
        # WARNING: The key signature here is critical! Used in ui.py to
        # separate out figure keywords and gridspec keywords.
        def _assign_scalar(key, value):
            if value is None:
                return
            if not np.isscalar(value):
                raise ValueError(f'Unexpected {key}={value!r}. Must be scalar.')
            value = units(value, 'em', 'in')
            setattr(self, '_' + key, value)
        hequal = _not_none(hequal, equal)
        wequal = _not_none(wequal, equal)
        _assign_scalar('left', left)
        _assign_scalar('right', right)
        _assign_scalar('bottom', bottom)
        _assign_scalar('top', top)
        _assign_scalar('hequal', hequal)
        _assign_scalar('wequal', wequal)
        _assign_scalar('panelpad', panelpad)
        _assign_scalar('outerpad', outerpad)
        _assign_scalar('innerpad', innerpad)

        # Assign vector args
        # NOTE: Here we employ obfuscation that skips 'panel' indices. So users could
        # still call self.update(wspace=[1, 2]) even if there is a right-axes panel
        # between each subplot. To control panel spaces users should instead pass
        # 'pad' or 'space' to panel_axes(), colorbar(), or legend() on creation.
        def _assign_vector(key, values, space):
            if values is None:
                return
            idxs = self._get_subplot_indices(key[0], space=space)
            nidxs = len(idxs)
            values = np.atleast_1d(values)
            if values.size == 1:
                values = np.repeat(values, nidxs)
            if values.size != nidxs:
                raise ValueError(f'Expected len({key}) == {nidxs}. Got {values.size}.')
            list_ = getattr(self, '_' + key)
            for i, value in enumerate(values):
                if value is None:
                    continue
                list_[idxs[i]] = value
        if pad is not None and not np.isscalar(pad):
            raise ValueError(f'Parameter pad={pad!r} must be scalar.')
        if space is not None and not np.isscalar(space):
            raise ValueError(f'Parameter space={space!r} must be scalar.')
        hpad = _not_none(hpad, pad)
        wpad = _not_none(wpad, pad)
        hpad = units(hpad, 'em', 'in')
        wpad = units(wpad, 'em', 'in')
        hspace = _not_none(hspace, space)
        wspace = _not_none(wspace, space)
        hspace = units(hspace, 'em', 'in')
        wspace = units(wspace, 'em', 'in')
        hratios = _not_none(hratios=hratios, height_ratios=height_ratios)
        wratios = _not_none(wratios=wratios, width_ratios=width_ratios)
        _assign_vector('hpad', hpad, space=True)
        _assign_vector('wpad', wpad, space=True)
        _assign_vector('hspace', hspace, space=True)
        _assign_vector('wspace', wspace, space=True)
        _assign_vector('hratios', hratios, space=False)
        _assign_vector('wratios', wratios, space=False)

    def get_geometry(self):
        """
        Return the total number of rows and columns in the grid.
        """
        return self.nrows, self.ncols

    def get_subplot_geometry(self):
        """
        Return the number of rows and columns allocated for "main" subplots
        and available with ``gridspec[...]`` indexing.
        """
        nrows, ncols = self.get_geometry()
        nrows_panels, ncols_panels = self.get_panel_geometry()
        return nrows - nrows_panels, ncols - ncols_panels

    def get_panel_geometry(self):
        """
        Return the number of rows and columns allocated for "panel" subplots
        and *not* available with ``gridspec[...]`` indexing.
        """
        nrows = sum(map(bool, self._hpanels))
        ncols = sum(map(bool, self._wpanels))
        return nrows, ncols

    def get_grid_positions(self, figure=None):
        """
        Return the subplot grid positions allowing for variable inter-subplot
        spacing and using physical units for the spacing terms.

        Note
        ----
        The physical units for positioning grid cells are converted from em-widths to
        inches when the `GridSpec` is instantiated. This means that subsequent changes
        to :rcraw:`font.size` will have no effect on the spaces. This is consistent
        with :rcraw:`font.size` having no effect on already-instantiated figures.
        """
        # Grab the figure size
        if not self.figure:
            self._figure = figure
        if not self.figure:
            raise RuntimeError('Figure must be assigned to gridspec.')
        if figure is not self.figure:
            raise RuntimeError('Cannot get positiona with non-gridspec figure.')
        fig = _not_none(figure, self.figure)
        figwidth, figheight = fig.get_size_inches()
        spacewidth, spaceheight = self.spacewidth, self.spaceheight
        panelwidth, panelheight = self.panelwidth, self.panelheight
        hratios, wratios = self.hratios, self.wratios
        hidxs = self._get_subplot_indices('h')
        widxs = self._get_subplot_indices('w')
        hsubplot = np.array([hratios[i] for i in hidxs])
        wsubplot = np.array([wratios[i] for i in widxs])

        # Scale the subplot slot ratios and keep the panel slots fixed
        hsubplot = (figheight - panelheight - spaceheight) * hsubplot / sum(hsubplot)
        wsubplot = (figwidth - panelwidth - spacewidth) * wsubplot / sum(wsubplot)
        hratios, hidxs = self.hratios, self._get_subplot_indices('h')
        for idx, ratio in zip(hidxs, hsubplot):
            hratios[idx] = ratio  # modify the main subplot ratios
        wratios, widxs = self.wratios, self._get_subplot_indices('w')
        for idx, ratio in zip(widxs, wsubplot):
            wratios[idx] = ratio

        # Calculate accumulated heights of columns
        norm = (figheight - spaceheight) / (figheight * sum(hratios))
        if norm < 0:
            raise RuntimeError(
                'Not enough room for axes. Try increasing the figure height or '
                "decreasing the 'top', 'bottom', or 'hspace' gridspec spaces."
            )
        cell_heights = [r * norm for r in hratios]
        sep_heights = [0] + [s / figheight for s in self.hspace]
        heights = np.cumsum(np.column_stack([sep_heights, cell_heights]).flat)

        # Calculate accumulated widths of rows
        norm = (figwidth - spacewidth) / (figwidth * sum(wratios))
        if norm < 0:
            raise RuntimeError(
                'Not enough room for axes. Try increasing the figure width or '
                "decreasing the 'left', 'right', or 'wspace' gridspec spaces."
            )
        cell_widths = [r * norm for r in wratios]
        sep_widths = [0] + [s / figwidth for s in self.wspace]
        widths = np.cumsum(np.column_stack([sep_widths, cell_widths]).flat)

        # Return the figure coordinates
        tops, bottoms = (1 - self.top / figheight - heights).reshape((-1, 2)).T
        lefts, rights = (self.left / figwidth + widths).reshape((-1, 2)).T
        return bottoms, tops, lefts, rights

    @docstring._snippet_manager
    def update(self, **kwargs):
        """
        Update the gridspec with arbitrary initialization keyword arguments
        and update the subplot positions.

        Parameters
        ----------
        %(gridspec.shared)s
        %(gridspec.vector)s
        %(gridspec.tight)s
        """
        # Apply positions to all axes
        # NOTE: This uses the current figure size to fix panel widths
        # and determine physical grid spacing.
        self._update_params(**kwargs)
        fig = self.figure
        for ax in fig.axes:
            if not isinstance(ax, maxes.SubplotBase):
                continue
            ss = ax.get_subplotspec().get_topmost_subplotspec()
            if ss.get_gridspec() is not self:  # should be impossible
                continue
            ax._reposition_subplot()
        fig.stale = True

    @property
    def figure(self):
        """
        The `proplot.figure.Figure` instance uniquely associated with this `GridSpec`.
        On assignment the gridspec parameters and figure size are updated.
        """
        return self._figure

    @figure.setter
    def figure(self, fig):
        from .figure import Figure
        if not isinstance(fig, Figure):
            raise ValueError('Figure must be a ProPlot figure.')
        self._figure = fig
        self._update_params(**fig._gridspec_params)
        fig._gridspec_params.clear()
        figsize = self._calc_figsize()
        if figsize is not None:
            fig.set_size_inches(figsize, internal=True, forward=False)
        else:
            self.update()

    # Delete attributes. Don't like having special setters and getters for some
    # settings and not others. Width and height ratios can be updated with update().
    # Also delete obsolete 'subplotpars' and built-in tight layout function.
    tight_layout = _disable_method('tight_layout')  # instead use custom tight layout
    subgridspec = _disable_method('subgridspec')  # instead use variable spaces
    get_width_ratios = _disable_method('get_width_ratios')
    get_height_ratios = _disable_method('get_height_ratios')
    set_width_ratios = _disable_method('set_width_ratios')
    set_height_ratios = _disable_method('set_height_ratios')
    get_subplot_params = _disable_method('get_subplot_params')
    locally_modified_subplot_params = _disable_method('locally_modified_subplot_params')

    # Make formerly public instance-level attributes immutable. Also redirect space
    # properties so they try to retrieve user settings then fallback to defaults.
    # NOTE: Do not document these since intended usage is internal and panel slot
    # obfuscation makes this confusing. For example gs.update(wspace=gs.wspace) in
    # presence of panels would yield error. For now the only supported introspection
    # is the __repr__. Probably no big deal... introspection not critical here.
    left = property(functools.partial(_get_current_space, key='left'), doc='')
    bottom = property(functools.partial(_get_current_space, key='bottom'), doc='')
    right = property(functools.partial(_get_current_space, key='right'), doc='')
    top = property(functools.partial(_get_current_space, key='top'), doc='')
    hspace = property(functools.partial(_get_current_space, key='hspace'), doc='')
    wspace = property(functools.partial(_get_current_space, key='wspace'), doc='')
    # Additional properties added for consistency
    nrows = property(lambda self: self._nrows, doc='')  # in case missing
    ncols = property(lambda self: self._ncols, doc='')  # ...
    hratios = property(lambda self: list(self._hratios))
    wratios = property(lambda self: list(self._wratios))
    hpad = property(lambda self: list(self._hpad))
    wpad = property(lambda self: list(self._wpad))
    # Hidden helper properties used to calculate figure size and subplot positions
    spaceheight = property(lambda self: self.bottom + self.top + sum(self.hspace))
    spacewidth = property(lambda self: self.left + self.right + sum(self.wspace))
    panelheight = property(
        lambda self: sum(r for i, r in enumerate(self.hratios) if self._hpanels[i])
    )
    panelwidth = property(
        lambda self: sum(r for i, r in enumerate(self.wratios) if self._wpanels[i])
    )
    subplotheight = property(
        lambda self: sum(r for i, r in enumerate(self.hratios) if not self._hpanels[i])
    )
    subplotwidth = property(
        lambda self: sum(r for i, r in enumerate(self.wratios) if not self._wpanels[i])
    )


class SubplotGrid(MutableSequence, list):
    """
    List-like object used to store subplots returned by
    `~proplot.figure.Figure.subplots`. 1D indexing uses the underlying list of
    `~proplot.axes.Axes` while 2D indexing uses the `~SubplotGrid.gridspec`.
    See `~SubplotGrid.__getitem__` for details.
    """
    def __repr__(self):
        if not self:
            return 'SubplotGrid(length=0)'
        length = len(self)
        nrows, ncols = self.gridspec.get_subplot_geometry()
        return f'SubplotGrid(nrows={nrows}, ncols={ncols}, length={length})'

    def __str__(self):
        return self.__repr__()

    def __len__(self):
        return list.__len__(self)

    def insert(self, key, value):  # required for MutableSequence
        value = self._validate_item(value, scalar=True)
        list.insert(self, key, value)

    def __init__(self, sequence=None, **kwargs):
        """
        Parameters
        ----------
        sequence : sequence
            A sequence of `proplot.axes.Axes` subplots or their children.

        See also
        --------
        proplot.figure.Figure.subplots
        proplot.ui.subplots
        """
        n = kwargs.pop('n', None)
        order = kwargs.pop('order', None)
        if n is not None or order is not None:
            warnings._warn_proplot(
                f'Ignoring n={n!r} and order={order!r}. As of v0.8 SubplotGrid '
                'handles 2D indexing by leveraging the subplotspec extents rather than '
                'directly emulating 2D array indexing. These arguments are no longer '
                'needed and will be removed in a future release.'
            )
        sequence = _not_none(sequence, [])
        sequence = self._validate_item(sequence, scalar=False)
        super().__init__(sequence, **kwargs)

    def __getattr__(self, attr):
        """
        Get a missing attribute. Simply redirects to the axes if the `SubplotGrid`
        is singleton and raises an error otherwise. This can be convenient for
        single-axes figures generated with `~proplot.figure.Figure.subplots`.
        """
        # Redirect to the axes
        if not self or attr[:1] == '_':
            return super().__getattribute__(attr)  # trigger default error
        if len(self) == 1:
            return getattr(self[0], attr)

        # Obscure deprecated behavior
        # WARNING: This is now deprecated! Instead we dynamically define a few
        # dedicated relevant commands that can be called from the grid (see below).
        warnings._warn_proplot(
            'Calling arbitrary axes methods from SubplotGrid was deprecated in v0.8 '
            'and will be removed in a future release. Please index the grid or loop '
            'over the grid instead.'
        )
        if not self:
            return None
        objs = tuple(getattr(ax, attr) for ax in self)  # may raise error
        if not any(map(callable, objs)):
            return objs[0] if len(self) == 1 else objs
        elif all(map(callable, objs)):
            @functools.wraps(objs[0])
            def _iterate_subplots(*args, **kwargs):
                result = []
                for func in objs:
                    result.append(func(*args, **kwargs))
                if len(self) == 1:
                    return result[0]
                elif all(res is None for res in result):
                    return None
                elif all(isinstance(res, paxes.Axes) for res in result):
                    return SubplotGrid(result, n=self._n, order=self._order)
                else:
                    return tuple(result)
            _iterate_subplots.__doc__ = inspect.getdoc(objs[0])
            return _iterate_subplots
        else:
            raise AttributeError(f'Found mixed types for attribute {attr!r}.')

    def __getitem__(self, key):
        """
        Get an axes.

        Parameters
        ----------
        key : int, slice, or 2-tuple
            The index. If 1D then the axes in the corresponding
            sublist are returned. If 2D then the axes that intersect
            the corresponding `~SubplotGrid.gridspec` slots are returned.

        Returns
        -------
        axs : `~proplot.axes.Axes` or `SubplotGrid`
            The axes. If the index included slices then
            another `SubplotGrid` is returned.

        Example
        -------
        >>> import proplot as pplt
        >>> fig, axs = pplt.subplots(nrows=3, ncols=3)
        >>> axs[3]  # the subplots in the second row, first column
        >>> axs[1, 2]  # the subplots in the second row, third column
        >>> axs[:, 0]  # a grid of subplots in the first column
        """
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]
        # List-style indexing
        if isinstance(key, (Integral, slice)):
            slices = isinstance(key, slice)
            objs = list.__getitem__(self, key)
        # Gridspec-style indexing
        elif (
            isinstance(key, tuple)
            and len(key) == 2
            and all(isinstance(ikey, (Integral, slice)) for ikey in key)
        ):
            # WARNING: Permit no-op slicing of empty grids here
            slices = any(isinstance(ikey, slice) for ikey in key)
            objs = []
            if self:
                gs = self.gridspec
                ss_key = gs._get_subplot_spec(key)  # obfuscates panels
                row1_key, col1_key = divmod(ss_key.num1, gs.ncols)
                row2_key, col2_key = divmod(ss_key.num2, gs.ncols)
            for ax in self:
                ss = ax._get_topmost_axes().get_subplotspec()
                row1, col1 = divmod(ss.num1, gs.ncols)
                row2, col2 = divmod(ss.num2, gs.ncols)
                inrow = row1_key <= row1 <= row2_key or row1_key <= row2 <= row2_key
                incol = col1_key <= col1 <= col2_key or col1_key <= col2 <= col2_key
                if inrow and incol:
                    objs.append(ax)
            if not slices and len(objs) == 1:
                objs = objs[0]
        else:
            raise IndexError(f'Invalid index {key!r}.')
        if isinstance(objs, list):
            return SubplotGrid(objs)
        else:
            return objs

    def __setitem__(self, key, value):
        """
        Add an axes.

        Parameters
        ----------
        key : int or slice
            The 1D index.
        value : `proplot.axes.Axes`
            The proplot subplot or its child or panel axes,
            or a sequence thereof if the index was a slice.
        """
        if isinstance(key, Integral):
            value = self._validate_item(value, scalar=True)
        elif isinstance(key, slice):
            value = self._validate_item(value, scalar=False)
        else:
            raise IndexError('Multi dimensional item assignment is not supported.')
        return super().__setitem__(key, value)  # could be list[:] = [1, 2, 3]

    def _validate_item(self, items, scalar=False):
        """
        Validate assignments. Accept diverse iterable inputs.
        """
        gridspec = None
        message = (
            'SubplotGrid can only be filled with ProPlot subplots '
            'belonging to the same GridSpec. Instead got {!r}.'
        )
        items = np.atleast_1d(items)
        if self:
            gridspec = self.gridspec  # compare against existing gridspec
        for item in items.flat:
            if not isinstance(item, paxes.Axes):
                raise ValueError(message.format(item))
            item = item._get_topmost_axes()
            if not isinstance(item, maxes.SubplotBase):
                raise ValueError(message.format(item))
            gs = item.get_subplotspec().get_gridspec()
            if not isinstance(gs, GridSpec) or (gridspec and gs is not gridspec):
                raise ValueError(message.format(gs))
            gridspec = gs
        if not scalar:
            items = tuple(items.flat)
        elif items.size == 1:
            items = items.flat[0]
        else:
            raise ValueError('Input must be a single ProPlot axes.')
        return items

    @docstring._snippet_manager
    def format(self, *args, **kwargs):
        """
        Call the ``format`` command for every axes in the grid.

        Parameters
        ----------
        %(axes.format)s
        **kwargs
            Passed to the projection-specific ``format`` command for each axes.
            Valid only if every axes in the grid belongs to the same class.

        Other parameters
        ----------------
        %(figure.format)s
        %(axes.rc)s

        See also
        --------
        proplot.axes.Axes.format
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        proplot.figure.Figure.format
        proplot.config.Configurator.context
        """
        for ax in self:
            ax.format(*args, **kwargs)

    @property
    def gridspec(self):
        """
        The `~proplot.gridspec.GridSpec` associated with the grid. This is used
        to resolve 2D indexing. See `~SubplotGrid.__getitem__` for details.
        """
        # Return the gridspec associatd with the grid
        if not self:
            raise ValueError('Unknown gridspec for empty SubplotGrid.')
        ax = self[0]
        ax = ax._get_topmost_axes()
        return ax.get_subplotspec().get_gridspec()

    @property
    def shape(self):
        """
        The shape of the `~proplot.gridspec.GridSpec` associated with the grid.
        See `~SubplotGrid.__getitem__` for details.
        """
        # NOTE: Considered deprecating this but on second thought since this is
        # a 2D array-like object it should definitely have a shape attribute.
        return self.gridspec.get_subplot_geometry()


def _add_grid_command(src, name):
    # Create the method
    def _grid_command(self, *args, **kwargs):
        objs = []
        for ax in self:
            obj = getattr(ax, name)(*args, **kwargs)
            objs.append(obj)
        return SubplotGrid(objs)

    # Clean the docstring
    cls = getattr(paxes, src)
    cmd = getattr(cls, name)
    doc = inspect.cleandoc(cmd.__doc__)  # dedents
    dot = doc.find('.')
    if dot != -1:
        doc = doc[:dot] + ' for every axes in the grid' + doc[dot:]
    doc = re.sub(
        r'^(Returns\n-------\n)(.+)(\n\s+)(.+)',
        r'\1SubplotGrid\2A grid of the resulting axes.',
        doc
    )

    # Apply the method
    _grid_command.__qualname__ = f'SubplotGrid.{name}'
    _grid_command.__name__ = name
    _grid_command.__doc__ = doc
    setattr(SubplotGrid, name, _grid_command)


# Dynamically add commands to generate twin or inset axes
# TODO: Add commands that plot the input data for every
# axes in the grid along a third dimension.
for _src, _name in (
    ('Axes', 'panel'),
    ('Axes', 'panel_axes'),
    ('Axes', 'inset'),
    ('Axes', 'inset_axes'),
    ('CartesianAxes', 'altx'),
    ('CartesianAxes', 'alty'),
    ('CartesianAxes', 'dualx'),
    ('CartesianAxes', 'dualy'),
    ('CartesianAxes', 'twinx'),
    ('CartesianAxes', 'twiny'),
):
    _add_grid_command(_src, _name)


# Deprecated
SubplotsContainer = warnings._rename_objs('0.8', SubplotsContainer=SubplotGrid)
