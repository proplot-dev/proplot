#!/usr/bin/env python3
"""
New gridspec and subplotspec classes.
"""
import matplotlib.axes as maxes
import matplotlib.gridspec as mgridspec
import numpy as np

from .config import rc
from .internals import ic  # noqa: F401
from .internals import _not_none, _version, _version_mpl
from .utils import units

__all__ = ['GridSpec', 'SubplotSpec']


def _default_space(key, share=0, pad=None):
    """
    Return suitable default spacing given a shared axes setting.
    """
    # Pull out sizes
    outerpad = _not_none(pad, rc['subplots.outerpad'])
    innerpad = _not_none(pad, rc['subplots.innerpad'])
    xtick = rc['xtick.major.size']
    ytick = rc['ytick.major.size']
    xtickpad = rc['xtick.major.pad']
    ytickpad = rc['ytick.major.pad']
    xticklabel = rc._scale_font(rc['xtick.labelsize'])
    yticklabel = 3 * rc._scale_font(rc['ytick.labelsize'])
    label = rc._scale_font(rc['axes.labelsize'])
    title = rc._scale_font(rc['axes.titlesize'])
    titlepad = rc['axes.titlepad']

    # Get suitable size for various spaces
    if key == 'left':
        space = units(outerpad) + (ytick + yticklabel + ytickpad + label) / 72
    elif key == 'right':
        space = units(outerpad)
    elif key == 'bottom':
        space = units(outerpad) + (xtick + xticklabel + xtickpad + label) / 72
    elif key == 'top':
        space = units(outerpad) + (title + titlepad) / 72
    elif key == 'wspace':
        space = units(innerpad) + ytick / 72
        if share < 3:
            space += (yticklabel + ytickpad) / 72
        if share < 1:
            space += label / 72
    elif key == 'hspace':
        space = units(innerpad) + (title + titlepad + xtick) / 72
        if share < 3:
            space += (xticklabel + xtickpad) / 72
        if share < 0:
            space += label / 72
    else:
        raise KeyError(f'Invalid space key {key!r}.')

    return space


def _calc_geometry(**kwargs):
    """
    Save arguments passed to `~proplot.ui.subplots`, calculates
    gridspec settings and figure size necessary for requested geometry, and
    returns keyword args necessary to reconstruct and modify this
    configuration. Note that `wspace`, `hspace`, `left`, `right`, `top`, and
    `bottom` always have fixed physical units, then we scale figure width,
    figure height, and width and height ratios to accommodate spaces.
    """
    # NOTE: In the future this will be replaced with a GeometryConfigurator
    # class that acts as the interface between the figure and the gridspec.
    # Dimensions and geometry
    refaspect = kwargs['refaspect']
    refxrange, refyrange = kwargs['refxrange'], kwargs['refyrange']
    refwidth, refheight = kwargs['refwidth'], kwargs['refheight']
    figwidth, figheight = kwargs['figwidth'], kwargs['figheight']
    nrows, ncols = kwargs['nrows'], kwargs['ncols']

    # Gridspec settings
    wspace, hspace = kwargs['wspace'], kwargs['hspace']
    wratios, hratios = kwargs['wratios'], kwargs['hratios']
    left, bottom = kwargs['left'], kwargs['bottom']
    right, top = kwargs['right'], kwargs['top']
    wequal, hequal, equal = kwargs['wequal'], kwargs['hequal'], kwargs['equal']

    # Panel string toggles, lists containing empty strings '' (indicating a
    # main axes), or one of 'l', 'r', 'b', 't' (indicating axes panels) or
    # 'f' (indicating figure panels)
    wpanels, hpanels = kwargs['wpanels'], kwargs['hpanels']

    # Checks, important now that we modify gridspec geometry
    if len(hratios) != nrows:
        raise ValueError(
            f'Expected {nrows} width ratios for {nrows} rows, '
            f'got {len(hratios)}.'
        )
    if len(wratios) != ncols:
        raise ValueError(
            f'Expected {ncols} width ratios for {ncols} columns, '
            f'got {len(wratios)}.'
        )
    if len(hspace) != nrows - 1:
        raise ValueError(
            f'Expected {nrows - 1} hspaces for {nrows} rows, '
            f'got {len(hspace)}.'
        )
    if len(wspace) != ncols - 1:
        raise ValueError(
            f'Expected {ncols - 1} wspaces for {ncols} columns, '
            f'got {len(wspace)}.'
        )
    if len(hpanels) != nrows:
        raise ValueError(
            f'Expected {nrows} hpanel toggles for {nrows} rows, '
            f'got {len(hpanels)}.'
        )
    if len(wpanels) != ncols:
        raise ValueError(
            f'Expected {ncols} wpanel toggles for {ncols} columns, '
            f'got {len(wpanels)}.'
        )

    # Get indices corresponding to main axes or main axes space slots
    idxs_ratios, idxs_space = [], []
    for panels in (hpanels, wpanels):
        # Ratio indices
        mask = np.array([bool(s) for s in panels])
        ratio_idxs, = np.where(~mask)
        idxs_ratios.append(ratio_idxs)
        # Space indices
        space_idxs = []
        for idx in ratio_idxs[:-1]:  # exclude last axes slot
            offset = 1
            while panels[idx + offset] not in 'rbf':  # main space next to this
                offset += 1
            space_idxs.append(idx + offset - 1)
        idxs_space.append(space_idxs)

    # Separate the panel and axes ratios
    hratios_main = [hratios[idx] for idx in idxs_ratios[0]]
    wratios_main = [wratios[idx] for idx in idxs_ratios[1]]
    hratios_panels = [r for idx, r in enumerate(hratios) if idx not in idxs_ratios[0]]
    wratios_panels = [r for idx, r in enumerate(wratios) if idx not in idxs_ratios[1]]
    hspace_main = [hspace[idx] for idx in idxs_space[0]]
    wspace_main = [wspace[idx] for idx in idxs_space[1]]

    # Reduced geometry
    nrows_main = len(hratios_main)
    ncols_main = len(wratios_main)

    # Get reference properties, account for panel slots in space and ratios
    # TODO: Shouldn't panel space be included in these calculations?
    (x1, x2), (y1, y2) = refxrange, refyrange
    dx, dy = x2 - x1 + 1, y2 - y1 + 1
    rwspace = sum(wspace_main[x1:x2])
    rhspace = sum(hspace_main[y1:y2])
    rwratio = ncols_main * sum(wratios_main[x1:x2 + 1]) / (dx * sum(wratios_main))
    rhratio = nrows_main * sum(hratios_main[y1:y2 + 1]) / (dy * sum(hratios_main))
    if rwratio == 0 or rhratio == 0:
        raise RuntimeError(
            f'Something went wrong, got wratio={rwratio!r} '
            f'and hratio={rhratio!r} for reference axes.'
        )
    if np.iterable(refaspect):
        refaspect = refaspect[0] / refaspect[1]

    # Determine figure and axes dims from input in width or height dimenion.
    # For e.g. common use case [[1,1,2,2],[0,3,3,0]], make sure we still scale
    # the reference axes like square even though takes two columns of gridspec!
    auto_width = figwidth is None and figheight is not None
    auto_height = figheight is None and figwidth is not None
    if figwidth is None and figheight is None:  # get stuff directly from axes
        if refwidth is None and refheight is None:
            refwidth = units(rc['subplots.refwidth'])
        if refheight is not None:
            auto_width = True
            refheight_all = (nrows_main * (refheight - rhspace)) / (dy * rhratio)
            figheight = refheight_all + top + bottom + sum(hspace) + sum(hratios_panels)
        if refwidth is not None:
            auto_height = True
            refwidth_all = (ncols_main * (refwidth - rwspace)) / (dx * rwratio)
            figwidth = refwidth_all + left + right + sum(wspace) + sum(wratios_panels)
        if refwidth is not None and refheight is not None:
            auto_width = auto_height = False
    else:
        if figheight is not None:
            refheight_all = figheight - top - bottom - sum(hspace) - sum(hratios_panels)
            refheight = (refheight_all * dy * rhratio) / nrows_main + rhspace
        if figwidth is not None:
            refwidth_all = figwidth - left - right - sum(wspace) - sum(wratios_panels)
            refwidth = (refwidth_all * dx * rwratio) / ncols_main + rwspace

    # Automatically figure dim that was not specified above
    if auto_height:
        refheight = refwidth / refaspect
        refheight_all = (nrows_main * (refheight - rhspace)) / (dy * rhratio)
        figheight = refheight_all + top + bottom + sum(hspace) + sum(hratios_panels)
    elif auto_width:
        refwidth = refheight * refaspect
        refwidth_all = (ncols_main * (refwidth - rwspace)) / (dx * rwratio)
        figwidth = refwidth_all + left + right + sum(wspace) + sum(wratios_panels)
    if refwidth_all < 0:
        raise ValueError(
            f'Not enough room for axes (would have width {refwidth_all}). '
            'Try using tight=False, increasing figure width, or decreasing '
            "'left', 'right', or 'wspace' spaces."
        )
    if refheight_all < 0:
        raise ValueError(
            f'Not enough room for axes (would have height {refheight_all}). '
            'Try using tight=False, increasing figure height, or decreasing '
            "'top', 'bottom', or 'hspace' spaces."
        )

    # Reconstruct the ratios array with physical units for subplot slots
    # The panel slots are unchanged because panels have fixed widths
    wratios_main = refwidth_all * np.array(wratios_main) / sum(wratios_main)
    hratios_main = refheight_all * np.array(hratios_main) / sum(hratios_main)
    for idx, ratio in zip(idxs_ratios[0], hratios_main):
        hratios[idx] = ratio
    for idx, ratio in zip(idxs_ratios[1], wratios_main):
        wratios[idx] = ratio

    # Convert margins to figure-relative coordinates
    left = left / figwidth
    bottom = bottom / figheight
    right = 1 - right / figwidth
    top = 1 - top / figheight

    # Constant spacing corrections
    if equal:  # do both
        wequal = hequal = True
    if wequal:
        wspace = [max(wspace) for _ in wspace]
    if hequal:
        hspace = [max(hspace) for _ in hspace]

    # Return gridspec keyword args
    gridspec_kw = {
        'ncols': ncols, 'nrows': nrows,
        'wspace': wspace, 'hspace': hspace,
        'width_ratios': wratios, 'height_ratios': hratios,
        'left': left, 'bottom': bottom, 'right': right, 'top': top,
    }

    figsize = (figwidth, figheight)
    return figsize, gridspec_kw, kwargs


class SubplotSpec(mgridspec.SubplotSpec):
    """
    Matplotlib `~matplotlib.gridspec.SubplotSpec` subclass that adds
    some helpful methods.
    """
    def __repr__(self):
        nrows, ncols, row1, row2, col1, col2 = self.get_rows_columns()
        return f'SubplotSpec({nrows}, {ncols}; {row1}:{row2}, {col1}:{col2})'

    def get_active_geometry(self):
        """
        Return the number of rows, number of columns, and 1d subplot
        location indices, ignoring rows and columns allocated for spaces.
        """
        nrows, ncols, row1, row2, col1, col2 = self.get_active_rows_columns()
        num1 = row1 * ncols + col1
        num2 = row2 * ncols + col2
        return nrows, ncols, num1, num2

    def get_active_rows_columns(self):
        """
        Return the number of rows, number of columns, first subplot row,
        last subplot row, first subplot column, and last subplot column,
        ignoring rows and columns allocated for spaces.
        """
        nrows, ncols = self.get_gridspec().get_geometry()
        row1, col1 = divmod(self.num1, ncols)
        if self.num2 is not None:
            row2, col2 = divmod(self.num2, ncols)
        else:
            row2 = row1
            col2 = col1
        return (
            (nrows + 1) // 2,
            (ncols + 1) // 2,
            row1 // 2,
            row2 // 2,
            col1 // 2,
            col2 // 2
        )

    def get_geometry(self):
        """
        Return the number of rows, number of columns, and 1d subplot
        location indices.
        """
        gridspec = self.get_gridspec()
        rows, cols = gridspec.get_geometry()
        return rows, cols, self.num1, self.num2

    def get_rows_columns(self):
        """
        Return the number of rows, number of columns, first subplot row,
        last subplot row, first subplot column, and last subplot column.
        """
        gridspec = self.get_gridspec()
        nrows, ncols = gridspec.get_geometry()
        row_start, col_start = divmod(self.num1, ncols)
        row_stop, col_stop = divmod(self.num2, ncols)
        return nrows, ncols, row_start, row_stop, col_start, col_stop

    @classmethod
    def _from_subplotspec(cls, subplotspec):
        """
        Translate a matplotlib subplotspec to proplot subplotspec.
        """
        return cls(subplotspec._gridspec, subplotspec.num1, subplotspec.num2)


class GridSpec(mgridspec.GridSpec):
    """
    Matplotlib `~matplotlib.gridspec.GridSpec` subclass that allows for grids
    with variable spacing between successive rows and columns of axes.
    This is done by drawing ``nrows * 2 + 1`` and ``ncols * 2 + 1``
    `~matplotlib.gridspec.GridSpec` rows and columns, setting `wspace`
    and `hspace` to ``0``, and masking out every other row and column
    of the `~matplotlib.gridspec.GridSpec`, so they act as "spaces".
    These "spaces" are then allowed to vary in width using the builtin
    `width_ratios` and `height_ratios` properties.

    Note
    ----
    In a future version, this class will natively support variable spacing
    between successive rows and columns without the obfuscation. It will also
    support specifying spaces in physical units via `~proplot.utils.units`.
    """
    def __repr__(self):  # do not show width and height ratios
        nrows, ncols = self.get_geometry()
        return f'GridSpec({nrows}, {ncols})'

    def __init__(self, figure, nrows=1, ncols=1, **kwargs):
        """
        Parameters
        ----------
        figure : `~proplot.figure.Figure`
            The figure instance filled by this gridspec. Unlike
            `~matplotlib.gridspec.GridSpec`, this argument is required.
        nrows, ncols : int, optional
            The number of rows and columns on the subplot grid.
        hspace, wspace : float or list of float
            The vertical and horizontal spacing between rows and columns of
            subplots, respectively. In `~proplot.ui.subplots`, ``wspace``
            and ``hspace`` are in physical units. When calling `GridSpec`
            directly, values are scaled relative to the average subplot
            height or width.

            If float, the spacing is identical between all rows and columns. If
            list of float, the length of the lists must equal ``nrows-1``
            and ``ncols-1``, respectively.
        height_ratios, width_ratios : list of float
            Ratios for the relative heights and widths for rows and columns
            of subplots, respectively. For example, ``width_ratios=(1,2)``
            scales a 2-column gridspec so that the second column is twice as
            wide as the first column.
        left, right, top, bottom : float or str
            Passed to `~matplotlib.gridspec.GridSpec`, denotes the margin
            positions in figure-relative coordinates.

        Other parameters
        ----------------
        **kwargs
            Passed to `~matplotlib.gridspec.GridSpec`.
        """
        self._nrows = nrows * 2 - 1  # used with get_geometry
        self._ncols = ncols * 2 - 1
        self._nrows_active = nrows
        self._ncols_active = ncols
        wratios, hratios, kwargs = self._spaces_as_ratios(**kwargs)
        super().__init__(
            self._nrows, self._ncols,
            hspace=0, wspace=0,  # replaced with "hidden" slots
            width_ratios=wratios, height_ratios=hratios,
            figure=figure, **kwargs
        )

    def __getitem__(self, key):
        """
        Magic obfuscation that renders `~matplotlib.gridspec.GridSpec`
        rows and columns designated as 'spaces' inaccessible.
        """
        nrows, ncols = self.get_geometry()
        nrows_active, ncols_active = self.get_active_geometry()
        if not isinstance(key, tuple):  # usage gridspec[1,2]
            num1, num2 = self._normalize(key, nrows_active * ncols_active)
        else:
            if len(key) == 2:
                k1, k2 = key
            else:
                raise ValueError(f'Invalid index {key!r}.')
            num1 = self._normalize(k1, nrows_active)
            num2 = self._normalize(k2, ncols_active)
            num1, num2 = np.ravel_multi_index((num1, num2), (nrows, ncols))
        num1 = self._positem(num1)
        num2 = self._positem(num2)
        return SubplotSpec(self, num1, num2)

    @staticmethod
    def _positem(size):
        """
        Account for negative indices.
        """
        if size < 0:
            # want -1 to stay -1, -2 becomes -3, etc.
            return 2 * (size + 1) - 1
        else:
            return size * 2

    @staticmethod
    def _normalize(key, size):
        """
        Transform gridspec index into standardized form.
        """
        if isinstance(key, slice):
            start, stop, _ = key.indices(size)
            if stop > start:
                return start, stop - 1
        else:
            if key < 0:
                key += size
            if 0 <= key < size:
                return key, key
        raise IndexError(f'Invalid index: {key} with size {size}.')

    def _spaces_as_ratios(
        self, hspace=None, wspace=None,  # spacing between axes
        height_ratios=None, width_ratios=None,
        **kwargs
    ):
        """
        For keyword argument usage, see `GridSpec`.
        """
        # Parse flexible input
        nrows, ncols = self.get_active_geometry()
        hratios = np.atleast_1d(_not_none(height_ratios, 1))
        wratios = np.atleast_1d(_not_none(width_ratios, 1))
        hspace = np.atleast_1d(_not_none(hspace, np.mean(hratios) * 0.10))  # relative
        wspace = np.atleast_1d(_not_none(wspace, np.mean(wratios) * 0.10))
        if len(wspace) == 1:
            wspace = np.repeat(wspace, (ncols - 1,))  # note: may be length 0
        if len(hspace) == 1:
            hspace = np.repeat(hspace, (nrows - 1,))
        if len(wratios) == 1:
            wratios = np.repeat(wratios, (ncols,))
        if len(hratios) == 1:
            hratios = np.repeat(hratios, (nrows,))

        # Verify input ratios and spacings
        # Translate height/width spacings, implement as extra columns/rows
        if len(hratios) != nrows:
            raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')
        if len(wratios) != ncols:
            raise ValueError(
                f'Got {ncols} columns, but {len(wratios)} wratios.'
            )
        if len(wspace) != ncols - 1:
            raise ValueError(
                f'Require {ncols-1} width spacings for {ncols} columns, '
                f'got {len(wspace)}.'
            )
        if len(hspace) != nrows - 1:
            raise ValueError(
                f'Require {nrows-1} height spacings for {nrows} rows, '
                f'got {len(hspace)}.'
            )

        # Assign spacing as ratios
        nrows, ncols = self.get_geometry()
        wratios_final = [None] * ncols
        wratios_final[::2] = list(wratios)
        if ncols > 1:
            wratios_final[1::2] = list(wspace)
        hratios_final = [None] * nrows
        hratios_final[::2] = list(hratios)
        if nrows > 1:
            hratios_final[1::2] = list(hspace)
        return wratios_final, hratios_final, kwargs  # bring extra kwargs back

    def get_margins(self):
        """
        Returns left, bottom, right, top values. Not sure why this method
        doesn't already exist on `~matplotlib.gridspec.GridSpec`.
        """
        return self.left, self.bottom, self.right, self.top

    def get_hspace(self):
        """
        Returns row ratios allocated for spaces.
        """
        return self.get_height_ratios()[1::2]

    def get_wspace(self):
        """
        Returns column ratios allocated for spaces.
        """
        return self.get_width_ratios()[1::2]

    def get_active_height_ratios(self):
        """
        Returns height ratios excluding slots allocated for spaces.
        """
        return self.get_height_ratios()[::2]

    def get_active_width_ratios(self):
        """
        Returns width ratios excluding slots allocated for spaces.
        """
        return self.get_width_ratios()[::2]

    def get_active_geometry(self):
        """
        Returns the number of active rows and columns, i.e. the rows and
        columns that aren't skipped by `~GridSpec.__getitem__`.
        """
        return self._nrows_active, self._ncols_active

    def update(self, **kwargs):
        """
        Update the gridspec with arbitrary initialization keyword arguments
        then *apply* those updates to every figure using this gridspec.
        The default `~matplotlib.gridspec.GridSpec.update` tries to update
        positions for axes on all active figures -- but this can fail after
        successive figure edits if it has been removed from the figure
        manager. ProPlot insists one gridspec per figure.

        Parameters
        ----------
        **kwargs
            Valid initialization keyword arguments. See `GridSpec`.
        """
        # Convert spaces to ratios
        wratios, hratios, kwargs = self._spaces_as_ratios(**kwargs)
        self.set_width_ratios(wratios)
        self.set_height_ratios(hratios)

        # Validate args
        nrows = kwargs.pop('nrows', None)
        ncols = kwargs.pop('ncols', None)
        nrows_current, ncols_current = self.get_active_geometry()
        if (
            nrows is not None and nrows != nrows_current
            or ncols is not None and ncols != ncols_current
        ):
            raise ValueError(
                f'Input geometry {(nrows, ncols)} does not match '
                f'current geometry {(nrows_current, ncols_current)}.'
            )
        self.left = kwargs.pop('left', None)
        self.right = kwargs.pop('right', None)
        self.bottom = kwargs.pop('bottom', None)
        self.top = kwargs.pop('top', None)
        if kwargs:
            raise ValueError(f'Unknown keyword arg(s): {kwargs}.')

        # Apply to figure and all axes
        fig = self.figure
        fig.subplotpars.update(self.left, self.bottom, self.right, self.top)
        for ax in fig.axes:
            if not isinstance(ax, maxes.SubplotBase):
                continue
            subplotspec = ax.get_subplotspec().get_topmost_subplotspec()
            if subplotspec.get_gridspec() is not self:
                continue
            if _version_mpl >= _version('3.4.0'):
                ax.set_position(ax.get_subplotspec().get_position(ax.figure))
            else:
                ax.update_params()
                ax.set_position(ax.figbox)  # equivalent to above
        fig.stale = True


class _GridSpecFromSubplotSpec(mgridspec.GridSpecFromSubplotSpec):
    """
    Subclass that generates `SubplotSpec` objects. Avoids a recent deprecation
    of `~matplotlib.gridspec.SubplotSpec.get_rows_columns`. This is currently only
    used to draw colorbars with partial span along subplot edges but will not be
    needed once proplot implements "edge stacks."
    """
    def __getitem__(self, key):
        subplotspec = super().__getitem__(key)
        return SubplotSpec._from_subplotspec(subplotspec)
