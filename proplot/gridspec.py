#!/usr/bin/env python3
"""
The gridspec and subplot grid classes used throughout proplot.
"""
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

__all__ = ["GridSpec", "SubplotGrid", "SubplotsContainer"]  # deprecated


# Gridspec vector arguments
# Valid for figure() and GridSpec()
_shared_docstring = """
left, right, top, bottom : unit-spec, default: None
    The fixed space between the subplots and the figure edge.
    %(units.em)s
    If ``None``, the space is determined automatically based on the tick and
    label settings. If :rcraw:`subplots.tight` is ``True`` or ``tight=True`` was
    passed to the figure, the space is determined by the tight layout algorithm.
"""
_scalar_docstring = """
wspace, hspace, space : unit-spec, default: None
    The fixed space between grid columns, rows, or both.
    %(units.em)s
    If ``None``, the space is determined automatically based on the font size and axis
    sharing settings. If :rcraw:`subplots.tight` is ``True`` or ``tight=True`` was
    passed to the figure, the space is determined by the tight layout algorithm.
"""
_vector_docstring = """
wspace, hspace, space : unit-spec or sequence, default: None
    The fixed space between grid columns, rows, and both, respectively. If
    float, string, or ``None``, this value is expanded into lists of length
    ``ncols - 1`` (for `wspace`) or length ``nrows - 1`` (for `hspace`). If
    a sequence, its length must match these lengths.
    %(units.em)s

    For elements equal to ``None``, the space is determined automatically based
    on the tick and label settings. If :rcraw:`subplots.tight` is ``True`` or
    ``tight=True`` was passed to the figure, the space is determined by the tight
    layout algorithm. For example, ``subplots(ncols=3, tight=True, wspace=(2, None))``
    fixes the space between columns 1 and 2 but lets the tight layout algorithm
    determine the space between columns 2 and 3.
wratios, hratios : float or sequence, optional
    Passed to `~proplot.gridspec.GridSpec`, denotes the width and height
    ratios for the subplot grid. Length of `wratios` must match the number
    of columns, and length of `hratios` must match the number of rows.
width_ratios, height_ratios
    Aliases for `wratios`, `hratios`. Included for
    consistency with `matplotlib.gridspec.GridSpec`.
wpad, hpad, pad : unit-spec or sequence, optional
    The tight layout padding between columns, rows, and both, respectively.
    Unlike ``space``, these control the padding between subplot content
    (including text, ticks, etc.) rather than subplot edges. As with
    ``space``, these can be scalars or arrays optionally containing ``None``.
    For elements equal to ``None``, the default is `innerpad`.
    %(units.em)s
"""
_tight_docstring = """
wequal, hequal, equal :  bool, default: :rc:`subplots.equalspace`
    Whether to make the tight layout algorithm apply equal spacing
    between columns, rows, or both.
wgroup, hgroup, group :  bool, default: :rc:`subplots.groupspace`
    Whether to make the tight layout algorithm just consider spaces between
    adjacent subplots instead of entire columns and rows of subplots.
outerpad : unit-spec, default: :rc:`subplots.outerpad`
    The scalar tight layout padding around the left, right, top, bottom figure edges.
    %(units.em)s
innerpad : unit-spec, default: :rc:`subplots.innerpad`
    The scalar tight layout padding between columns and rows. Synonymous with `pad`.
    %(units.em)s
panelpad : unit-spec, default: :rc:`subplots.panelpad`
    The scalar tight layout padding between subplots and their panels,
    colorbars, and legends and between "stacks" of these objects.
    %(units.em)s
"""
docstring._snippet_manager["gridspec.shared"] = _shared_docstring
docstring._snippet_manager["gridspec.scalar"] = _scalar_docstring
docstring._snippet_manager["gridspec.vector"] = _vector_docstring
docstring._snippet_manager["gridspec.tight"] = _tight_docstring


def _disable_method(attr):
    """
    Disable the inherited method.
    """

    def _dummy_method(*args):
        raise RuntimeError(f"Method {attr}() is disabled on proplot gridspecs.")

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
            nrows, ncols, num1, num2 = self._get_geometry()
        except (IndexError, ValueError, AttributeError):
            return "SubplotSpec(unknown)"
        else:
            return f"SubplotSpec(nrows={nrows}, ncols={ncols}, index=({num1}, {num2}))"

    def _get_geometry(self):
        """
        Return the geometry and scalar indices relative to the "unhidden" non-panel
        geometry. May trigger error if this is in a "hidden" panel slot.
        """
        gs = self.get_gridspec()
        num1, num2 = self.num1, self.num2
        if isinstance(gs, GridSpec):
            nrows, ncols = gs.get_geometry()
            num1, num2 = gs._decode_indices(num1, num2)  # may trigger error
        return nrows, ncols, num1, num2

    def _get_rows_columns(self, ncols=None):
        """
        Return the row and column indices. The resulting indices include
        "hidden" panel rows and columns. See `GridSpec.get_grid_positions`.
        """
        # NOTE: Sort of confusing that this doesn't have 'total' in name but that
        # is by analogy with get_grid_positions(). This is used for grid positioning.
        gs = self.get_gridspec()
        if isinstance(gs, GridSpec):
            ncols = _not_none(ncols, gs.ncols_total)
        else:
            ncols = _not_none(ncols, gs.ncols)
        row1, col1 = divmod(self.num1, ncols)
        row2, col2 = divmod(self.num2, ncols)
        return row1, row2, col1, col2

    def get_position(self, figure, return_all=False):
        # Silent override. Older matplotlib versions can create subplots
        # with negative heights and widths that crash on instantiation.
        # Instead better to dynamically adjust the bounding box and hope
        # that subsequent adjustments will correct the subplot position.
        gs = self.get_gridspec()
        if isinstance(gs, GridSpec):
            nrows, ncols = gs.get_total_geometry()
        else:
            nrows, ncols = gs.get_geometry()
        rows, cols = np.unravel_index([self.num1, self.num2], (nrows, ncols))
        bottoms, tops, lefts, rights = gs.get_grid_positions(figure)
        bottom = bottoms[rows].min()
        top = max(bottom, tops[rows].max())
        left = lefts[cols].min()
        right = max(left, rights[cols].max())
        bbox = mtransforms.Bbox.from_extents(left, bottom, right, top)
        if return_all:
            return bbox, rows[0], cols[0], nrows, ncols
        else:
            return bbox


class GridSpec(mgridspec.GridSpec):
    """
    A `~matplotlib.gridspec.GridSpec` subclass that permits variable spacing
    between successive rows and columns and hides "panel slots" from indexing.
    """

    def __repr__(self):
        nrows, ncols = self.get_geometry()
        prows, pcols = self.get_panel_geometry()
        params = {"nrows": nrows, "ncols": ncols}
        if prows:
            params["nrows_panel"] = prows
        if pcols:
            params["ncols_panel"] = pcols
        params = ", ".join(f"{key}={value!r}" for key, value in params.items())
        return f"GridSpec({params})"

    def __getattr__(self, attr):
        # Redirect to private 'layout' attributes that are fragile w.r.t.
        # matplotlib version. Cannot set these by calling super().__init__()
        # because we make spacing arguments non-settable properties.
        if "layout" in attr:
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

        See also
        --------
        proplot.ui.figure
        proplot.figure.Figure
        proplot.ui.subplots
        proplot.figure.Figure.subplots
        proplot.figure.Figure.add_subplots
        matplotlib.gridspec.GridSpec

        Important
        ---------
        Adding axes panels, axes or figure colorbars, and axes or figure legends
        quietly augments the gridspec geometry by inserting "panel slots". However,
        subsequently indexing the gridspec with ``gs[num]`` or ``gs[row, col]`` will
        ignore the "panel slots". This permits adding new subplots by passing
        ``gs[num]`` or ``gs[row, col]`` to `~proplot.figure.Figure.add_subplot`
        even in the presence of panels (see `~GridSpec.__getitem__` for details).
        This also means that each `GridSpec` is `~proplot.figure.Figure`-specific,
        i.e. it can only be used once (if you are working with `GridSpec` instances
        manually and want the same geometry for multiple figures, you must create
        a copy with `GridSpec.copy` before working on the subsequent figure).
        """
        # Fundamental GridSpec properties
        self._nrows_total = nrows
        self._ncols_total = ncols
        self._left = None
        self._right = None
        self._bottom = None
        self._top = None
        self._hspace_total = [None] * (nrows - 1)
        self._wspace_total = [None] * (ncols - 1)
        self._hratios_total = [1] * nrows
        self._wratios_total = [1] * ncols
        self._left_default = None
        self._right_default = None
        self._bottom_default = None
        self._top_default = None
        self._hspace_total_default = [None] * (nrows - 1)
        self._wspace_total_default = [None] * (ncols - 1)
        self._figure = None  # initial state

        # Capture rc settings used for default spacing
        # NOTE: This is consistent with conversion of 'em' units to inches on gridspec
        # instantiation. In general it seems strange for future changes to rc settings
        # to magically update an existing gridspec layout. This also may improve draw
        # time as manual or auto figure resizes repeatedly call get_grid_positions().
        scales = {"in": 0, "inout": 0.5, "out": 1, None: 1}
        self._xtickspace = scales[rc["xtick.direction"]] * rc["xtick.major.size"]
        self._ytickspace = scales[rc["ytick.direction"]] * rc["ytick.major.size"]
        self._xticklabelspace = (
            _fontsize_to_pt(rc["xtick.labelsize"]) + rc["xtick.major.pad"]
        )  # noqa: E501
        self._yticklabelspace = (
            2 * _fontsize_to_pt(rc["ytick.labelsize"]) + rc["ytick.major.pad"]
        )  # noqa: E501
        self._labelspace = _fontsize_to_pt(rc["axes.labelsize"]) + rc["axes.labelpad"]
        self._titlespace = _fontsize_to_pt(rc["axes.titlesize"]) + rc["axes.titlepad"]

        # Tight layout and panel-related properties
        # NOTE: The wpanels and hpanels contain empty strings '' (indicating main axes),
        # or one of 'l', 'r', 'b', 't' (indicating axes panels) or 'f' (figure panels)
        outerpad = _not_none(kwargs.pop("outerpad", None), rc["subplots.outerpad"])
        innerpad = _not_none(kwargs.pop("innerpad", None), rc["subplots.innerpad"])
        panelpad = _not_none(kwargs.pop("panelpad", None), rc["subplots.panelpad"])
        pad = _not_none(kwargs.pop("pad", None), innerpad)  # alias of innerpad
        self._outerpad = units(outerpad, "em", "in")
        self._innerpad = units(innerpad, "em", "in")
        self._panelpad = units(panelpad, "em", "in")
        self._hpad_total = [units(pad, "em", "in")] * (nrows - 1)
        self._wpad_total = [units(pad, "em", "in")] * (ncols - 1)
        self._hequal = rc["subplots.equalspace"]
        self._wequal = rc["subplots.equalspace"]
        self._hgroup = rc["subplots.groupspace"]
        self._wgroup = rc["subplots.groupspace"]
        self._hpanels = [""] * nrows  # axes and figure panel identification
        self._wpanels = [""] * ncols
        self._fpanels = {  # array representation of figure panel spans
            "left": np.empty((0, nrows), dtype=bool),
            "right": np.empty((0, nrows), dtype=bool),
            "bottom": np.empty((0, ncols), dtype=bool),
            "top": np.empty((0, ncols), dtype=bool),
        }
        self._update_params(pad=pad, **kwargs)

    def __getitem__(self, key):
        """
        Get a `~matplotlib.gridspec.SubplotSpec`. "Hidden" slots allocated for axes
        panels, colorbars, and legends are ignored. For example, given a gridspec with
        2 subplot rows, 3 subplot columns, and a "panel" row between the subplot rows,
        calling ``gs[1, 1]`` returns a `~matplotlib.gridspec.SubplotSpec` corresponding
        to the central subplot on the second row rather than a "panel" slot.
        """
        return self._make_subplot_spec(key, includepanels=False)

    def _make_subplot_spec(self, key, includepanels=False):
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
            extra = "for gridspec" if axis is None else f"along axis {axis}"
            raise IndexError(f"Invalid index {key} {extra} with size {size}.")

        # Normalize the indices
        if includepanels:
            nrows, ncols = self.get_total_geometry()
        else:
            nrows, ncols = self.get_geometry()
        if not isinstance(key, tuple):  # usage gridspec[1,2]
            num1, num2 = _normalize_index(key, nrows * ncols)
        elif len(key) == 2:
            k1, k2 = key
            num1 = _normalize_index(k1, nrows, axis=0)
            num2 = _normalize_index(k2, ncols, axis=1)
            num1, num2 = np.ravel_multi_index((num1, num2), (nrows, ncols))
        else:
            raise ValueError(f"Invalid index {key!r}.")

        # Return the subplotspec
        if not includepanels:
            num1, num2 = self._encode_indices(num1, num2)
        return _SubplotSpec(self, num1, num2)

    def _encode_indices(self, *args, which=None):
        """
        Convert indices from the "unhidden" gridspec geometry into indices for the
        total geometry. If `which` is not passed these should be flattened indices.
        """
        nums = []
        idxs = self._get_indices(which)
        for arg in args:
            try:
                nums.append(idxs[arg])
            except (IndexError, TypeError):
                raise ValueError(f"Invalid gridspec index {arg}.")
        return nums[0] if len(nums) == 1 else nums

    def _decode_indices(self, *args, which=None):
        """
        Convert indices from the total geometry into the "unhidden" gridspec
        geometry. If `which` is not passed these should be flattened indices.
        """
        nums = []
        idxs = self._get_indices(which)
        for arg in args:
            try:
                nums.append(idxs.index(arg))
            except ValueError:
                raise ValueError(f"Invalid gridspec index {arg}.")
        return nums[0] if len(nums) == 1 else nums

    def _filter_indices(self, key, panel=False):
        """
        Filter the vector attribute for "unhidden" or "hidden" slots.
        """
        # NOTE: Currently this is just used for unused internal properties,
        # defined for consistency with the properties ending in "total".
        # These may be made public in a future version.
        which = key[0]
        space = "space" in key or "pad" in key
        idxs = self._get_indices(which=which, space=space, panel=panel)
        vector = getattr(self, key + "_total")
        return [vector[i] for i in idxs]

    def _get_indices(self, which=None, space=False, panel=False):
        """
        Get the indices associated with "unhidden" or "hidden" slots.
        """
        if which:
            panels = getattr(self, f"_{which}panels")
        else:
            panels = [h + w for h, w in itertools.product(self._hpanels, self._wpanels)]
        if not space:
            idxs = [i for i, p in enumerate(panels) if p]
        else:
            idxs = [
                i
                for i, (p1, p2) in enumerate(zip(panels[:-1], panels[1:]))
                if p1 == p2 == "f"
                or p1 in ("l", "t")
                and p2 in ("l", "t", "")
                or p1 in ("r", "b", "")
                and p2 in ("r", "b")
            ]
        if not panel:
            length = len(panels) - 1 if space else len(panels)
            idxs = [i for i in range(length) if i not in idxs]
        return idxs

    def _modify_subplot_geometry(self, newrow=None, newcol=None):
        """
        Update the axes subplot specs by inserting rows and columns as specified.
        """
        fig = self.figure
        ncols = self._ncols_total - int(newcol is not None)  # previous columns
        inserts = (newrow, newrow, newcol, newcol)
        for ax in fig._iter_axes(hidden=True, children=True):
            # Get old index
            # NOTE: Endpoints are inclusive, not exclusive!
            if not isinstance(ax, maxes.SubplotBase):
                continue
            gs = ax.get_subplotspec().get_gridspec()
            ss = ax.get_subplotspec().get_topmost_subplotspec()
            # Get a new subplotspec
            coords = list(ss._get_rows_columns(ncols=ncols))
            for i in range(4):
                if inserts[i] is not None and coords[i] >= inserts[i]:
                    coords[i] += 1
            row1, row2, col1, col2 = coords
            key1 = slice(row1, row2 + 1)
            key2 = slice(col1, col2 + 1)
            ss_new = self._make_subplot_spec((key1, key2), includepanels=True)
            # Apply new subplotspec
            # NOTE: We should only have one possible level of GridSpecFromSubplotSpec
            # nesting -- from making side colorbars with length less than 1.
            if ss is ax.get_subplotspec():
                ax.set_subplotspec(ss_new)
            elif ss is getattr(gs, "_subplot_spec", None):
                gs._subplot_spec = ss_new
            else:
                raise RuntimeError("Unexpected GridSpecFromSubplotSpec nesting.")
            ax._reposition_subplot()

    def _parse_panel_arg(self, side, arg):
        """
        Return the indices associated with a new figure panel on the specified side.
        Try to find room in the current mosaic of figure panels.
        """
        # Add a subplot panel. Index depends on the side
        # NOTE: This always "stacks" new panels on old panels
        if isinstance(arg, maxes.SubplotBase) and isinstance(arg, paxes.Axes):
            slot = side[0]
            ss = arg.get_subplotspec().get_topmost_subplotspec()
            offset = len(arg._panel_dict[side]) + 1
            row1, row2, col1, col2 = ss._get_rows_columns()
            if side in ("left", "right"):
                iratio = col1 - offset if side == "left" else col2 + offset
                start, stop = row1, row2
            else:
                iratio = row1 - offset if side == "top" else row2 + offset
                start, stop = col1, col2

        # Add a figure panel. Index depends on the side and the input 'span'
        # NOTE: Here the 'span' indices start at '1' by analogy with add_subplot()
        # integers and with main subplot numbers. Also *ignores panel slots*.
        # NOTE: This only "stacks" panels if requested slots are filled. Slots are
        # tracked with figure panel array (a boolean mask where each row corresponds
        # to a panel, moving toward the outside, and True indicates a slot is filled).
        elif (
            arg is None
            or isinstance(arg, Integral)
            or np.iterable(arg)
            and all(isinstance(_, Integral) for _ in arg)
        ):
            slot = "f"
            array = self._fpanels[side]
            nacross = (
                self._ncols_total if side in ("left", "right") else self._nrows_total
            )  # noqa: E501
            npanels, nalong = array.shape
            arg = np.atleast_1d(_not_none(arg, (1, nalong)))
            if arg.size not in (1, 2):
                raise ValueError(
                    f"Invalid span={arg!r}. Must be scalar or 2-tuple of coordinates."
                )  # noqa: E501
            if any(s < 1 or s > nalong for s in arg):
                raise ValueError(
                    f"Invalid span={arg!r}. Coordinates must satisfy 1 <= c <= {nalong}."
                )  # noqa: E501
            start, stop = arg[0] - 1, arg[-1]  # non-inclusive starting at zero
            iratio = -1 if side in ("left", "top") else nacross  # default values
            for i in range(npanels):  # possibly use existing panel slot
                if not any(array[i, start:stop]):
                    array[i, start:stop] = True
                    if side in ("left", "top"):  # descending moves us closer to 0
                        iratio = npanels - 1 - i  # index in ratios array
                    else:  # descending array moves us closer to nacross - 1
                        iratio = nacross - (npanels - i)  # index in ratios array
                    break
            if iratio == -1 or iratio == nacross:  # no slots so we must add to array
                iarray = np.zeros((1, nalong), dtype=bool)
                iarray[0, start:stop] = True
                array = np.concatenate((array, iarray), axis=0)
                self._fpanels[side] = array  # replace array
            which = "h" if side in ("left", "right") else "w"
            start, stop = self._encode_indices(start, stop - 1, which=which)

        else:
            raise ValueError(f"Invalid panel argument {arg!r}.")

        # Return subplotspec indices
        # NOTE: Convert using the lengthwise indices
        return slot, iratio, slice(start, stop + 1)

    def _insert_panel_slot(
        self,
        side,
        arg,
        *,
        share=None,
        width=None,
        space=None,
        pad=None,
        filled=False,
    ):
        """
        Insert a panel slot into the existing gridspec. The `side` is the panel side
        and the `arg` is either an axes instance or the figure row-column span.
        """
        # Parse input args and get user-input properties, default properties
        fig = self.figure
        if fig is None:
            raise RuntimeError("Figure must be assigned to gridspec.")
        if side not in ("left", "right", "bottom", "top"):
            raise ValueError(f"Invalid side {side}.")
        slot, idx, span = self._parse_panel_arg(side, arg)
        pad = units(pad, "em", "in")
        space = units(space, "em", "in")
        width = units(width, "in")
        share = False if filled else share if share is not None else True
        which = "w" if side in ("left", "right") else "h"
        panels = getattr(self, f"_{which}panels")
        pads = getattr(self, f"_{which}pad_total")  # no copies!
        ratios = getattr(self, f"_{which}ratios_total")
        spaces = getattr(self, f"_{which}space_total")
        spaces_default = getattr(self, f"_{which}space_total_default")
        new_outer_slot = idx in (-1, len(panels))
        new_inner_slot = not new_outer_slot and panels[idx] != slot

        # Retrieve default spaces
        # NOTE: Cannot use 'wspace' and 'hspace' for top and right colorbars because
        # that adds an unnecessary tick space. So bypass _get_default_space totally.
        pad_default = (
            self._panelpad
            if slot != "f"
            or side in ("left", "top")
            and panels[0] == "f"
            or side in ("right", "bottom")
            and panels[-1] == "f"
            else self._innerpad
        )
        inner_space_default = (
            _not_none(pad, pad_default)
            if side in ("top", "right")
            else self._get_default_space(
                "hspace_total" if side == "bottom" else "wspace_total",
                title=False,  # no title between subplot and panel
                share=3 if share else 0,  # space for main subplot labels
                pad=_not_none(pad, pad_default),
            )
        )
        outer_space_default = self._get_default_space(
            (
                "bottom"
                if not share and side == "top"
                else "left" if not share and side == "right" else side
            ),
            title=True,  # room for titles deflected above panels
            pad=self._outerpad if new_outer_slot else self._innerpad,
        )
        if new_inner_slot:
            outer_space_default += self._get_default_space(
                "hspace_total" if side in ("bottom", "top") else "wspace_total",
                share=None,  # use external share setting
                pad=0,  # use no additional padding
            )
        width_default = units(
            rc["colorbar.width" if filled else "subplots.panelwidth"], "in"
        )

        # Adjust space, ratio, and panel indicator arrays
        # If slot exists, overwrite width, pad, space if they were provided by the user
        # If slot does not exist, modify gemoetry and add insert new spaces
        attr = "ncols" if side in ("left", "right") else "nrows"
        idx_offset = int(side in ("top", "left"))
        idx_inner_space = idx - int(side in ("bottom", "right"))  # inner colorbar space
        idx_outer_space = idx - int(side in ("top", "left"))  # outer colorbar space
        if new_outer_slot or new_inner_slot:
            idx += idx_offset
            idx_inner_space += idx_offset
            idx_outer_space += idx_offset
            newcol, newrow = (idx, None) if attr == "ncols" else (None, idx)
            setattr(self, f"_{attr}_total", 1 + getattr(self, f"_{attr}_total"))
            panels.insert(idx, slot)
            ratios.insert(idx, _not_none(width, width_default))
            pads.insert(idx_inner_space, _not_none(pad, pad_default))
            spaces.insert(idx_inner_space, space)
            spaces_default.insert(idx_inner_space, inner_space_default)
            if new_inner_slot:
                spaces_default.insert(idx_outer_space, outer_space_default)
            else:
                setattr(self, f"_{side}_default", outer_space_default)
        else:
            newrow = newcol = None
            spaces_default[idx_inner_space] = inner_space_default
            if width is not None:
                ratios[idx] = width
            if pad is not None:
                pads[idx_inner_space] = pad
            if space is not None:
                spaces[idx_inner_space] = space

        # Update the figure and axes and return a SubplotSpec
        # NOTE: For figure panels indices are determined by user-input spans.
        self._modify_subplot_geometry(newrow, newcol)
        figsize = self._update_figsize()
        if figsize is not None:
            fig.set_size_inches(figsize, internal=True, forward=False)
        else:
            self.update()
        key = (span, idx) if side in ("left", "right") else (idx, span)
        ss = self._make_subplot_spec(key, includepanels=True)  # bypass obfuscation
        return ss, share

    def _get_space(self, key):
        """
        Return the currently active vector inner space or scalar outer space
        accounting for both default values and explicit user overrides.
        """
        # NOTE: Default panel spaces should have been filled by _insert_panel_slot.
        # They use 'panelpad' and the panel-local 'share' setting. This function
        # instead fills spaces between subplots depending on sharing setting.
        fig = self.figure
        if not fig:
            raise ValueError("Figure must be assigned to get grid positions.")
        attr = f"_{key}"  # user-specified
        attr_default = f"_{key}_default"  # default values
        value = getattr(self, attr)
        value_default = getattr(self, attr_default)
        if key in ("left", "right", "bottom", "top"):
            if value_default is None:
                value_default = self._get_default_space(key)
                setattr(self, attr_default, value_default)
            return _not_none(value, value_default)
        elif key in ("wspace_total", "hspace_total"):
            result = []
            for i, (val, val_default) in enumerate(zip(value, value_default)):
                if val_default is None:
                    val_default = self._get_default_space(key)
                    value_default[i] = val_default
                result.append(_not_none(val, val_default))
            return result
        else:
            raise ValueError(f"Unknown space parameter {key!r}.")

    def _get_default_space(self, key, pad=None, share=None, title=True):
        """
        Return suitable default scalar inner or outer space given a shared axes
        setting. This is only relevant when "tight layout" is disabled.
        """
        # NOTE: Internal spacing args are stored in inches to simplify the
        # get_grid_positions() calculations.
        fig = self.figure
        if fig is None:
            raise RuntimeError("Figure must be assigned.")
        if key == "right":
            pad = _not_none(pad, self._outerpad)
            space = 0
        elif key == "top":
            pad = _not_none(pad, self._outerpad)
            space = self._titlespace if title else 0
        elif key == "left":
            pad = _not_none(pad, self._outerpad)
            space = self._labelspace + self._yticklabelspace + self._ytickspace
        elif key == "bottom":
            pad = _not_none(pad, self._outerpad)
            space = self._labelspace + self._xticklabelspace + self._xtickspace
        elif key == "wspace_total":
            pad = _not_none(pad, self._innerpad)
            share = _not_none(share, fig._sharey, 0)
            space = self._ytickspace
            if share < 3:
                space += self._yticklabelspace
            if share < 1:
                space += self._labelspace
        elif key == "hspace_total":
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
            raise ValueError(f"Invalid space key {key!r}.")
        return pad + space / 72

    def _get_tight_space(self, w):
        """
        Get tight layout spaces between the input subplot rows or columns.
        """
        # Get constants
        fig = self.figure
        if not fig:
            return
        if w == "w":
            x, y = "xy"
            group = self._wgroup
            nacross = self.nrows_total
            space = self.wspace_total
            pad = self.wpad_total
        else:
            x, y = "yx"
            group = self._hgroup
            nacross = self.ncols_total
            space = self.hspace_total
            pad = self.hpad_total

        # Iterate along each row or column space
        axs = tuple(fig._iter_axes(hidden=True, children=False))
        space = list(space)  # a copy
        ralong = np.array([ax._range_subplotspec(x) for ax in axs])
        racross = np.array([ax._range_subplotspec(y) for ax in axs])
        for i, (s, p) in enumerate(zip(space, pad)):
            # Find axes that abutt aginst this row or column space
            groups = []
            for j in range(nacross):  # e.g. each row
                # Get the indices for axes that meet this row or column edge.
                # NOTE: Rigorously account for empty and overlapping slots here
                filt = (racross[:, 0] <= j) & (j <= racross[:, 1])
                if sum(filt) < 2:
                    continue  # no interface
                ii = i
                idx1 = idx2 = np.array(())
                while ii >= 0 and idx1.size == 0:
                    filt1 = ralong[:, 1] == ii  # i.e. r / b edge abutts against this
                    (idx1,) = np.where(filt & filt1)
                    ii -= 1
                ii = i + 1
                while ii <= len(space) and idx2.size == 0:
                    filt2 = ralong[:, 0] == ii  # i.e. l / t edge abutts against this
                    (idx2,) = np.where(filt & filt2)
                    ii += 1
                # Put axes into unique groups and store as (l, r) or (b, t) pairs.
                axs1, axs2 = [axs[_] for _ in idx1], [axs[_] for _ in idx2]
                if x != "x":  # order bottom-to-top
                    axs1, axs2 = axs2, axs1
                for group1, group2 in groups:
                    if any(_ in group1 for _ in axs1) or any(_ in group2 for _ in axs2):
                        group1.update(axs1)
                        group2.update(axs2)
                        break
                else:
                    if axs1 and axs2:
                        groups.append((set(axs1), set(axs2)))  # form new group
            # Determing the spaces using cached tight bounding boxes
            # NOTE: Set gridspec space to zero if there are no adjacent edges
            if not group:
                groups = [
                    (
                        set(ax for (group1, _) in groups for ax in group1),
                        set(ax for (_, group2) in groups for ax in group2),
                    )
                ]
            margins = []
            for group1, group2 in groups:
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
        if ratio == "auto":
            return
        elif ratio == "equal":
            ratio = 1
        elif isinstance(ratio, str):
            raise RuntimeError(f"Unknown aspect ratio mode {ratio!r}.")
        else:
            ratio = 1 / ratio

        # Compare to current aspect after scaling by data ratio
        # Noat matplotlib 3.2.0 expanded get_data_ratio to work for all axis scales:
        # https://github.com/matplotlib/matplotlib/commit/87c742b99dc6b9a190f8c89bc6256ced72f5ab80  # noqa: E501
        aspect = ratio / ax.get_data_ratio()
        if fig._refaspect is not None:
            return  # fixed by user
        if np.isclose(aspect, fig._refaspect_default):
            return  # close enough to the default aspect
        fig._refaspect_default = aspect

        # Update the layout
        figsize = self._update_figsize()
        if not fig._is_same_size(figsize):
            fig.set_size_inches(figsize, internal=True)

    def _auto_layout_tight(self, renderer):
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
        hspace = self._get_tight_space("h")
        wspace = self._get_tight_space("w")
        if self._hequal:
            idxs = self._get_indices("h", space=True)
            space = max(hspace[i] for i in idxs)
            for i in idxs:
                hspace[i] = space
        if self._wequal:
            idxs = self._get_indices("w", space=True)
            space = max(wspace[i] for i in idxs)
            for i in idxs:
                wspace[i] = space
        self._hspace_total_default = hspace
        self._wspace_total_default = wspace

        # Update the layout
        # NOTE: fig.set_size_inches() always updates the gridspec to enforce fixed
        # spaces (necessary since native position coordinates are figure-relative)
        # and to enforce fixed panel ratios. So only self.update() if we skip resize.
        figsize = self._update_figsize()
        if not fig._is_same_size(figsize):
            fig.set_size_inches(figsize, internal=True)
        else:
            self.update()

    def _update_figsize(self):
        """
        Return an updated auto layout figure size accounting for the
        gridspec and figure parameters. May or may not need to be applied.
        """
        fig = self.figure
        if fig is None:  # drawing before subplots are added?
            return
        ax = fig._subplot_dict.get(fig._refnum, None)
        if ax is None:  # drawing before subplots are added?
            return
        ss = ax.get_subplotspec().get_topmost_subplotspec()
        y1, y2, x1, x2 = ss._get_rows_columns()
        refhspace = sum(self.hspace_total[y1:y2])
        refwspace = sum(self.wspace_total[x1:x2])
        refhpanel = sum(
            self.hratios_total[i] for i in range(y1, y2 + 1) if self._hpanels[i]
        )  # noqa: E501
        refwpanel = sum(
            self.wratios_total[i] for i in range(x1, x2 + 1) if self._wpanels[i]
        )  # noqa: E501
        refhsubplot = sum(
            self.hratios_total[i] for i in range(y1, y2 + 1) if not self._hpanels[i]
        )  # noqa: E501
        refwsubplot = sum(
            self.wratios_total[i] for i in range(x1, x2 + 1) if not self._wpanels[i]
        )  # noqa: E501

        # Get the reference sizes
        # NOTE: The sizing arguments should have been normalized already
        figwidth, figheight = fig._figwidth, fig._figheight
        refwidth, refheight = fig._refwidth, fig._refheight
        refaspect = _not_none(fig._refaspect, fig._refaspect_default)
        if refheight is None and figheight is None:
            if figwidth is not None:
                gridwidth = figwidth - self.spacewidth - self.panelwidth
                refwidth = gridwidth * refwsubplot / self.gridwidth
            if refwidth is not None:  # WARNING: do not change to elif!
                refheight = refwidth / refaspect
            else:
                raise RuntimeError("Figure size arguments are all missing.")
        if refwidth is None and figwidth is None:
            if figheight is not None:
                gridheight = figheight - self.spaceheight - self.panelheight
                refheight = gridheight * refhsubplot / self.gridheight
            if refheight is not None:
                refwidth = refheight * refaspect
            else:
                raise RuntimeError("Figure size arguments are all missing.")

        # Get the auto figure size. Might trigger 'not enough room' error later
        # NOTE: For e.g. [[1, 1, 2, 2], [0, 3, 3, 0]] we make sure to still scale the
        # reference axes like a square even though takes two columns of gridspec.
        if refheight is not None:
            refheight -= refhspace + refhpanel
            gridheight = refheight * self.gridheight / refhsubplot
            figheight = gridheight + self.spaceheight + self.panelheight
        if refwidth is not None:
            refwidth -= refwspace + refwpanel
            gridwidth = refwidth * self.gridwidth / refwsubplot
            figwidth = gridwidth + self.spacewidth + self.panelwidth

        # Return the figure size
        figsize = (figwidth, figheight)
        if all(np.isfinite(figsize)):
            return figsize
        else:
            warnings._warn_proplot(f"Auto resize failed. Invalid figsize {figsize}.")

    def _update_params(
        self,
        *,
        left=None,
        bottom=None,
        right=None,
        top=None,
        wspace=None,
        hspace=None,
        space=None,
        wpad=None,
        hpad=None,
        pad=None,
        wequal=None,
        hequal=None,
        equal=None,
        wgroup=None,
        hgroup=None,
        group=None,
        outerpad=None,
        innerpad=None,
        panelpad=None,
        hratios=None,
        wratios=None,
        width_ratios=None,
        height_ratios=None,
    ):
        """
        Update the user-specified properties.
        """

        # Assign scalar args
        # WARNING: The key signature here is critical! Used in ui.py to
        # separate out figure keywords and gridspec keywords.
        def _assign_scalar(key, value, convert=True):
            if value is None:
                return
            if not np.isscalar(value):
                raise ValueError(f"Unexpected {key}={value!r}. Must be scalar.")
            if convert:
                value = units(value, "em", "in")
            setattr(self, f"_{key}", value)

        hequal = _not_none(hequal, equal)
        wequal = _not_none(wequal, equal)
        hgroup = _not_none(hgroup, group)
        wgroup = _not_none(wgroup, group)
        _assign_scalar("left", left)
        _assign_scalar("right", right)
        _assign_scalar("bottom", bottom)
        _assign_scalar("top", top)
        _assign_scalar("panelpad", panelpad)
        _assign_scalar("outerpad", outerpad)
        _assign_scalar("innerpad", innerpad)
        _assign_scalar("hequal", hequal, convert=False)
        _assign_scalar("wequal", wequal, convert=False)
        _assign_scalar("hgroup", hgroup, convert=False)
        _assign_scalar("wgroup", wgroup, convert=False)

        # Assign vector args
        # NOTE: Here we employ obfuscation that skips 'panel' indices. So users could
        # still call self.update(wspace=[1, 2]) even if there is a right-axes panel
        # between each subplot. To control panel spaces users should instead pass
        # 'pad' or 'space' to panel_axes(), colorbar(), or legend() on creation.
        def _assign_vector(key, values, space):
            if values is None:
                return
            idxs = self._get_indices(key[0], space=space)
            nidxs = len(idxs)
            values = np.atleast_1d(values)
            if values.size == 1:
                values = np.repeat(values, nidxs)
            if values.size != nidxs:
                raise ValueError(f"Expected len({key}) == {nidxs}. Got {values.size}.")
            list_ = getattr(self, f"_{key}_total")
            for i, value in enumerate(values):
                if value is None:
                    continue
                list_[idxs[i]] = value

        if pad is not None and not np.isscalar(pad):
            raise ValueError(f"Parameter pad={pad!r} must be scalar.")
        if space is not None and not np.isscalar(space):
            raise ValueError(f"Parameter space={space!r} must be scalar.")
        hpad = _not_none(hpad, pad)
        wpad = _not_none(wpad, pad)
        hpad = units(hpad, "em", "in")
        wpad = units(wpad, "em", "in")
        hspace = _not_none(hspace, space)
        wspace = _not_none(wspace, space)
        hspace = units(hspace, "em", "in")
        wspace = units(wspace, "em", "in")
        hratios = _not_none(hratios=hratios, height_ratios=height_ratios)
        wratios = _not_none(wratios=wratios, width_ratios=width_ratios)
        _assign_vector("hpad", hpad, space=True)
        _assign_vector("wpad", wpad, space=True)
        _assign_vector("hspace", hspace, space=True)
        _assign_vector("wspace", wspace, space=True)
        _assign_vector("hratios", hratios, space=False)
        _assign_vector("wratios", wratios, space=False)

    @docstring._snippet_manager
    def copy(self, **kwargs):
        """
        Return a copy of the `GridSpec` with the `~proplot.figure.Figure`-specific
        "panel slots" removed. This can be useful if you want to draw multiple
        figures with the same geometry. Properties are inherited from this
        `GridSpec` by default but can be changed by passing keyword arguments.

        Parameters
        ----------
        %(gridspec.shared)s
        %(gridspec.vector)s
        %(gridspec.tight)s

        See also
        --------
        GridSpec.update
        """
        # WARNING: For some reason copy.copy() fails. Updating e.g. wpanels
        # and hpanels on the copy also updates this object. No idea why.
        nrows, ncols = self.get_geometry()
        gs = GridSpec(nrows, ncols)
        hidxs = self._get_indices("h")
        widxs = self._get_indices("w")
        gs._hratios_total = [self._hratios_total[i] for i in hidxs]
        gs._wratios_total = [self._wratios_total[i] for i in widxs]
        hidxs = self._get_indices("h", space=True)
        widxs = self._get_indices("w", space=True)
        gs._hpad_total = [self._hpad_total[i] for i in hidxs]
        gs._wpad_total = [self._wpad_total[i] for i in widxs]
        gs._hspace_total = [self._hspace_total[i] for i in hidxs]
        gs._wspace_total = [self._wspace_total[i] for i in widxs]
        gs._hspace_total_default = [self._hspace_total_default[i] for i in hidxs]
        gs._wspace_total_default = [self._wspace_total_default[i] for i in widxs]
        for key in (
            "left",
            "right",
            "bottom",
            "top",
            "labelspace",
            "titlespace",
            "xtickspace",
            "ytickspace",
            "xticklabelspace",
            "yticklabelspace",
            "outerpad",
            "innerpad",
            "panelpad",
            "hequal",
            "wequal",
        ):
            value = getattr(self, "_" + key)
            setattr(gs, "_" + key, value)
        gs.update(**kwargs)
        return gs

    def get_geometry(self):
        """
        Return the number of "unhidden" non-panel rows and columns in the grid
        (see `GridSpec` for details).

        See also
        --------
        GridSpec.get_panel_geometry
        GridSpec.get_total_geometry
        """
        nrows, ncols = self.get_total_geometry()
        nrows_panels, ncols_panels = self.get_panel_geometry()
        return nrows - nrows_panels, ncols - ncols_panels

    def get_panel_geometry(self):
        """
        Return the number of "hidden" panel rows and columns in the grid
        (see `GridSpec` for details).

        See also
        --------
        GridSpec.get_geometry
        GridSpec.get_total_geometry
        """
        nrows = sum(map(bool, self._hpanels))
        ncols = sum(map(bool, self._wpanels))
        return nrows, ncols

    def get_total_geometry(self):
        """
        Return the total number of "unhidden" and "hidden" rows and columns
        in the grid (see `GridSpec` for details).

        See also
        --------
        GridSpec.get_geometry
        GridSpec.get_panel_geometry
        GridSpec.get_grid_positions
        """
        return self._nrows_total, self._ncols_total

    def get_grid_positions(self, figure=None):
        """
        Return the subplot grid positions allowing for variable inter-subplot
        spacing and using physical units for the spacing terms. The resulting
        positions include "hidden" panel rows and columns.

        Note
        ----
        The physical units for positioning grid cells are converted from em-widths to
        inches when the `GridSpec` is instantiated. This means that subsequent changes
        to :rcraw:`font.size` will have no effect on the spaces. This is consistent
        with :rcraw:`font.size` having no effect on already-instantiated figures.

        See also
        --------
        GridSpec.get_total_geometry
        """
        # Grab the figure size
        if not self.figure:
            self._figure = figure
        if not self.figure:
            raise RuntimeError("Figure must be assigned to gridspec.")
        if figure is not self.figure:
            raise RuntimeError(
                f"Input figure {figure} does not match gridspec figure {self.figure}."
            )  # noqa: E501
        fig = _not_none(figure, self.figure)
        figwidth, figheight = fig.get_size_inches()
        spacewidth, spaceheight = self.spacewidth, self.spaceheight
        panelwidth, panelheight = self.panelwidth, self.panelheight
        hratios, wratios = self.hratios_total, self.wratios_total
        hidxs, widxs = self._get_indices("h"), self._get_indices("w")

        # Scale the subplot slot ratios and keep the panel slots fixed
        hsubplot = np.array([hratios[i] for i in hidxs])
        wsubplot = np.array([wratios[i] for i in widxs])
        hsubplot = (figheight - panelheight - spaceheight) * hsubplot / np.sum(hsubplot)
        wsubplot = (figwidth - panelwidth - spacewidth) * wsubplot / np.sum(wsubplot)
        for idx, ratio in zip(hidxs, hsubplot):
            hratios[idx] = ratio  # modify the main subplot ratios
        for idx, ratio in zip(widxs, wsubplot):
            wratios[idx] = ratio

        # Calculate accumulated heights of columns
        norm = (figheight - spaceheight) / (figheight * sum(hratios))
        if norm < 0:
            raise RuntimeError(
                "Not enough room for axes. Try increasing the figure height or "
                "decreasing the 'top', 'bottom', or 'hspace' gridspec spaces."
            )
        cell_heights = [r * norm for r in hratios]
        sep_heights = [0] + [s / figheight for s in self.hspace_total]
        heights = np.cumsum(np.column_stack([sep_heights, cell_heights]).flat)

        # Calculate accumulated widths of rows
        norm = (figwidth - spacewidth) / (figwidth * sum(wratios))
        if norm < 0:
            raise RuntimeError(
                "Not enough room for axes. Try increasing the figure width or "
                "decreasing the 'left', 'right', or 'wspace' gridspec spaces."
            )
        cell_widths = [r * norm for r in wratios]
        sep_widths = [0] + [s / figwidth for s in self.wspace_total]
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

        See also
        --------
        GridSpec.copy
        """
        # Apply positions to all axes
        # NOTE: This uses the current figure size to fix panel widths
        # and determine physical grid spacing.
        self._update_params(**kwargs)
        fig = self.figure
        if fig is None:
            return
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
        The `proplot.figure.Figure` uniquely associated with this `GridSpec`.
        On assignment the gridspec parameters and figure size are updated.

        See also
        --------
        proplot.gridspec.SubplotGrid.figure
        proplot.figure.Figure.gridspec
        """
        return self._figure

    @figure.setter
    def figure(self, fig):
        from .figure import Figure

        if not isinstance(fig, Figure):
            raise ValueError("Figure must be a proplot figure.")
        if self._figure and self._figure is not fig:
            raise ValueError(
                "Cannot use the same gridspec for multiple figures. "
                "Please use gridspec.copy() to make a copy."
            )
        self._figure = fig
        self._update_params(**fig._gridspec_params)
        fig._gridspec_params.clear()
        figsize = self._update_figsize()
        if figsize is not None:
            fig.set_size_inches(figsize, internal=True, forward=False)
        else:
            self.update()

    # Delete attributes. Don't like having special setters and getters for some
    # settings and not others. Width and height ratios can be updated with update().
    # Also delete obsolete 'subplotpars' and built-in tight layout function.
    tight_layout = _disable_method("tight_layout")  # instead use custom tight layout
    subgridspec = _disable_method("subgridspec")  # instead use variable spaces
    get_width_ratios = _disable_method("get_width_ratios")
    get_height_ratios = _disable_method("get_height_ratios")
    set_width_ratios = _disable_method("set_width_ratios")
    set_height_ratios = _disable_method("set_height_ratios")
    get_subplot_params = _disable_method("get_subplot_params")
    locally_modified_subplot_params = _disable_method("locally_modified_subplot_params")

    # Immutable helper properties used to calculate figure size and subplot positions
    # NOTE: The spaces are auto-filled with defaults wherever user left them unset
    gridheight = property(lambda self: sum(self.hratios))
    gridwidth = property(lambda self: sum(self.wratios))
    panelheight = property(lambda self: sum(self.hratios_panel))
    panelwidth = property(lambda self: sum(self.wratios_panel))
    spaceheight = property(lambda self: self.bottom + self.top + sum(self.hspace_total))
    spacewidth = property(lambda self: self.left + self.right + sum(self.wspace_total))

    # Geometry properties. These are included for consistency with get_geometry
    # functions (would be really confusing if self.nrows, self.ncols disagree).
    nrows = property(
        lambda self: self._nrows_total - sum(map(bool, self._hpanels)), doc=""
    )  # noqa: E501
    ncols = property(
        lambda self: self._ncols_total - sum(map(bool, self._wpanels)), doc=""
    )  # noqa: E501
    nrows_panel = property(lambda self: sum(map(bool, self._hpanels)))
    ncols_panel = property(lambda self: sum(map(bool, self._wpanels)))
    nrows_total = property(lambda self: self._nrows_total)
    ncols_total = property(lambda self: self._ncols_total)

    # Make formerly public instance-level attributes immutable and redirect space
    # properties so they try to retrieve user settings then fallback to defaults.
    # NOTE: These are undocumented for the time being. Generally properties should
    # be changed with update() and introspection not really necessary.
    left = property(lambda self: self._get_space("left"))
    bottom = property(lambda self: self._get_space("bottom"))
    right = property(lambda self: self._get_space("right"))
    top = property(lambda self: self._get_space("top"))
    hratios = property(lambda self: self._filter_indices("hratios", panel=False))
    wratios = property(lambda self: self._filter_indices("wratios", panel=False))
    hratios_panel = property(lambda self: self._filter_indices("hratios", panel=True))
    wratios_panel = property(lambda self: self._filter_indices("wratios", panel=True))
    hratios_total = property(lambda self: list(self._hratios_total))
    wratios_total = property(lambda self: list(self._wratios_total))
    hspace = property(lambda self: self._filter_indices("hspace", panel=False))
    wspace = property(lambda self: self._filter_indices("wspace", panel=False))
    hspace_panel = property(lambda self: self._filter_indices("hspace", panel=True))
    wspace_panel = property(lambda self: self._filter_indices("wspace", panel=True))
    hspace_total = property(lambda self: self._get_space("hspace_total"))
    wspace_total = property(lambda self: self._get_space("wspace_total"))
    hpad = property(lambda self: self._filter_indices("hpad", panel=False))
    wpad = property(lambda self: self._filter_indices("wpad", panel=False))
    hpad_panel = property(lambda self: self._filter_indices("hpad", panel=True))
    wpad_panel = property(lambda self: self._filter_indices("wpad", panel=True))
    hpad_total = property(lambda self: list(self._hpad_total))
    wpad_total = property(lambda self: list(self._wpad_total))


class SubplotGrid(MutableSequence, list):
    """
    List-like, array-like object used to store subplots returned by
    `~proplot.figure.Figure.subplots`. 1D indexing uses the underlying list of
    `~proplot.axes.Axes` while 2D indexing uses the `~SubplotGrid.gridspec`.
    See `~SubplotGrid.__getitem__` for details.
    """

    def __repr__(self):
        if not self:
            return "SubplotGrid(length=0)"
        length = len(self)
        nrows, ncols = self.gridspec.get_geometry()
        return f"SubplotGrid(nrows={nrows}, ncols={ncols}, length={length})"

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
        proplot.ui.subplots
        proplot.figure.Figure.subplots
        proplot.figure.Figure.add_subplots
        """
        n = kwargs.pop("n", None)
        order = kwargs.pop("order", None)
        if n is not None or order is not None:
            warnings._warn_proplot(
                f"Ignoring n={n!r} and order={order!r}. As of v0.8 SubplotGrid "
                "handles 2D indexing by leveraging the subplotspec extents rather than "
                "directly emulating 2D array indexing. These arguments are no longer "
                "needed and will be removed in a future release."
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
        if not self or attr[:1] == "_":
            return super().__getattribute__(attr)  # trigger default error
        if len(self) == 1:
            return getattr(self[0], attr)

        # Obscure deprecated behavior
        # WARNING: This is now deprecated! Instead we dynamically define a few
        # dedicated relevant commands that can be called from the grid (see below).
        import functools

        warnings._warn_proplot(
            "Calling arbitrary axes methods from SubplotGrid was deprecated in v0.8 "
            "and will be removed in a future release. Please index the grid or loop "
            "over the grid instead."
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
            raise AttributeError(f"Found mixed types for attribute {attr!r}.")

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
        axs : proplot.axes.Axes or SubplotGrid
            The axes. If the index included slices then
            another `SubplotGrid` is returned.

        Example
        -------
        >>> import proplot as pplt
        >>> fig, axs = pplt.subplots(nrows=3, ncols=3)
        >>> axs[5]  # the subplot in the second row, third column
        >>> axs[1, 2]  # the subplot in the second row, third column
        >>> axs[:, 0]  # a SubplotGrid containing the subplots in the first column
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
                ss_key = gs._make_subplot_spec(key)  # obfuscates panels
                row1_key, col1_key = divmod(ss_key.num1, gs.ncols)
                row2_key, col2_key = divmod(ss_key.num2, gs.ncols)
            for ax in self:
                ss = ax._get_topmost_axes().get_subplotspec().get_topmost_subplotspec()
                row1, col1 = divmod(ss.num1, gs.ncols)
                row2, col2 = divmod(ss.num2, gs.ncols)
                inrow = row1_key <= row1 <= row2_key or row1_key <= row2 <= row2_key
                incol = col1_key <= col1 <= col2_key or col1_key <= col2 <= col2_key
                if inrow and incol:
                    objs.append(ax)
            if not slices and len(objs) == 1:  # accounts for overlapping subplots
                objs = objs[0]
        else:
            raise IndexError(f"Invalid index {key!r}.")
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
            raise IndexError("Multi dimensional item assignment is not supported.")
        return super().__setitem__(key, value)  # could be list[:] = [1, 2, 3]

    @classmethod
    def _add_command(cls, src, name):
        """
        Add a `SubplotGrid` method that iterates through axes methods.
        """

        # Create the method
        def _grid_command(self, *args, **kwargs):
            objs = []
            for ax in self:
                obj = getattr(ax, name)(*args, **kwargs)
                objs.append(obj)
            return SubplotGrid(objs)

        # Clean the docstring
        cmd = getattr(src, name)
        doc = inspect.cleandoc(cmd.__doc__)  # dedents
        dot = doc.find(".")
        if dot != -1:
            doc = doc[:dot] + " for every axes in the grid" + doc[dot:]
        doc = re.sub(
            r"^(Returns\n-------\n)(.+)(\n\s+)(.+)",
            r"\1SubplotGrid\2A grid of the resulting axes.",
            doc,
        )

        # Apply the method
        _grid_command.__qualname__ = f"SubplotGrid.{name}"
        _grid_command.__name__ = name
        _grid_command.__doc__ = doc
        setattr(cls, name, _grid_command)

    def _validate_item(self, items, scalar=False):
        """
        Validate assignments. Accept diverse iterable inputs.
        """
        gridspec = None
        message = (
            "SubplotGrid can only be filled with proplot subplots "
            "belonging to the same GridSpec. Instead got {}."
        )
        items = np.atleast_1d(items)
        if self:
            gridspec = self.gridspec  # compare against existing gridspec
        for item in items.flat:
            if not isinstance(item, paxes.Axes):
                raise ValueError(message.format(f"the object {item!r}"))
            item = item._get_topmost_axes()
            if not isinstance(item, maxes.SubplotBase):
                raise ValueError(message.format(f"the axes {item!r}"))
            gs = item.get_subplotspec().get_topmost_subplotspec().get_gridspec()
            if not isinstance(gs, GridSpec):
                raise ValueError(message.format(f"the GridSpec {gs!r}"))
            if gridspec and gs is not gridspec:
                raise ValueError(message.format("at least two different GridSpecs"))
            gridspec = gs
        if not scalar:
            items = tuple(items.flat)
        elif items.size == 1:
            items = items.flat[0]
        else:
            raise ValueError("Input must be a single proplot axes.")
        return items

    @docstring._snippet_manager
    def format(self, **kwargs):
        """
        Call the ``format`` command for the `~SubplotGrid.figure`
        and every axes in the grid.

        Parameters
        ----------
        %(axes.format)s
        **kwargs
            Passed to the projection-specific ``format`` command for each axes.
            Valid only if every axes in the grid belongs to the same class.

        Other parameters
        ----------------
        %(figure.format)s
        %(cartesian.format)s
        %(polar.format)s
        %(geo.format)s
        %(rc.format)s

        See also
        --------
        proplot.axes.Axes.format
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        proplot.figure.Figure.format
        proplot.config.Configurator.context
        """
        self.figure.format(axs=self, **kwargs)

    @property
    def figure(self):
        """
        The `proplot.figure.Figure` uniquely associated with this `SubplotGrid`.
        This is used with the `SubplotGrid.format` command.

        See also
        --------
        proplot.gridspec.GridSpec.figure
        proplot.gridspec.SubplotGrid.gridspec
        proplot.figure.Figure.subplotgrid
        """
        return self.gridspec.figure

    @property
    def gridspec(self):
        """
        The `~proplot.gridspec.GridSpec` uniquely associated with this `SubplotGrid`.
        This is used to resolve 2D indexing. See `~SubplotGrid.__getitem__` for details.

        See also
        --------
        proplot.figure.Figure.gridspec
        proplot.gridspec.SubplotGrid.figure
        proplot.gridspec.SubplotGrid.shape
        """
        # Return the gridspec associatd with the grid
        if not self:
            raise ValueError("Unknown gridspec for empty SubplotGrid.")
        ax = self[0]
        ax = ax._get_topmost_axes()
        return ax.get_subplotspec().get_topmost_subplotspec().get_gridspec()

    @property
    def shape(self):
        """
        The shape of the `~proplot.gridspec.GridSpec` associated with the grid.
        See `~SubplotGrid.__getitem__` for details.

        See also
        --------
        proplot.gridspec.SubplotGrid.gridspec
        """
        # NOTE: Considered deprecating this but on second thought since this is
        # a 2D array-like object it should definitely have a shape attribute.
        return self.gridspec.get_geometry()


# Dynamically add commands to generate twin or inset axes
# TODO: Add commands that plot the input data for every
# axes in the grid along a third dimension.
for _src, _name in (
    (paxes.Axes, "panel"),
    (paxes.Axes, "panel_axes"),
    (paxes.Axes, "inset"),
    (paxes.Axes, "inset_axes"),
    (paxes.CartesianAxes, "altx"),
    (paxes.CartesianAxes, "alty"),
    (paxes.CartesianAxes, "dualx"),
    (paxes.CartesianAxes, "dualy"),
    (paxes.CartesianAxes, "twinx"),
    (paxes.CartesianAxes, "twiny"),
):
    SubplotGrid._add_command(_src, _name)

# Deprecated
SubplotsContainer = warnings._rename_objs("0.8.0", SubplotsContainer=SubplotGrid)
