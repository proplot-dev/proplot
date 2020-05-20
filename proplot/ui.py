#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures.
"""
import numpy as np
import functools
import inspect
import matplotlib.pyplot as plt
from . import constructor
from . import axes as paxes
from . import figure as pfigure
from . import gridspec as pgridspec
from .config import rc
from .utils import units
from .internals import ic  # noqa: F401
from .internals import warnings, _not_none

__all__ = [
    'close', 'show', 'subplots', 'SubplotsContainer', 'subplot_grid',
]

# Width or (width, height) dimensions for common journal specifications
JOURNAL_SPECS = {
    'aaas1': '5.5cm',
    'aaas2': '12cm',
    'agu1': ('95mm', '115mm'),
    'agu2': ('190mm', '115mm'),
    'agu3': ('95mm', '230mm'),
    'agu4': ('190mm', '230mm'),
    'ams1': 3.2,
    'ams2': 4.5,
    'ams3': 5.5,
    'ams4': 6.5,
    'nat1': '89mm',
    'nat2': '183mm',
    'pnas1': '8.7cm',
    'pnas2': '11.4cm',
    'pnas3': '17.8cm',
}


def close(*args, **kwargs):
    """
    Pass the input arguments to `matplotlib.pyplot.close`. This is included
    so you don't have to import `~matplotlib.pyplot`.
    """
    plt.close(*args, **kwargs)


def show():
    """
    Call `matplotlib.pyplot.show`. This is included so you don't have to import
    `~matplotlib.pyplot`. Note this command should *not be necessary* if you
    are working in an iPython session and :rcraw:`matplotlib` is non-empty --
    when you create a new figure, it will be automatically displayed.
    """
    plt.show()


def _journals(journal):
    """
    Return the width and height corresponding to the given journal.
    """
    # Get dimensions for figure from common journals.
    value = JOURNAL_SPECS.get(journal, None)
    if value is None:
        raise ValueError(
            f'Unknown journal figure size specifier {journal!r}. '
            'Current options are: '
            + ', '.join(map(repr, JOURNAL_SPECS.keys()))
        )
    # Return width, and optionally also the height
    width, height = None, None
    try:
        width, height = value
    except (TypeError, ValueError):
        width = value
    return width, height


def _axes_dict(naxs, value, kw=False, default=None):
    """
    Return a dictionary that looks like ``{1: value1, 2: value2, ...}`` or
    ``{1: {key1: value1, ...}, 2: {key2: value2, ...}, ...}`` for storing
    standardized axes-specific properties or keyword args.
    """
    # First build up dictionary
    # 1) 'string' or {1:'string1', (2,3):'string2'}
    if not kw:
        if np.iterable(value) and not isinstance(value, (str, dict)):
            value = {num + 1: item for num, item in enumerate(value)}
        elif not isinstance(value, dict):
            value = {range(1, naxs + 1): value}
    # 2) {'prop':value} or {1:{'prop':value1}, (2,3):{'prop':value2}}
    else:
        nested = [isinstance(value, dict) for value in value.values()]
        if not any(nested):  # any([]) == False
            value = {range(1, naxs + 1): value.copy()}
        elif not all(nested):
            raise ValueError(
                'Pass either of dictionary of key value pairs or '
                'a dictionary of dictionaries of key value pairs.'
            )
    # Then *unfurl* keys that contain multiple axes numbers, i.e. are meant
    # to indicate properties for multiple axes at once
    kwargs = {}
    for nums, item in value.items():
        nums = np.atleast_1d(nums)
        for num in nums.flat:
            if not kw:
                kwargs[num] = item
            else:
                kwargs[num] = item.copy()
    # Fill with default values
    for num in range(1, naxs + 1):
        if num not in kwargs:
            if kw:
                kwargs[num] = {}
            else:
                kwargs[num] = default
    # Verify numbers
    if {*range(1, naxs + 1)} != {*kwargs.keys()}:
        raise ValueError(
            f'Have {naxs} axes, but {value!r} has properties for axes '
            + ', '.join(map(repr, sorted(kwargs))) + '.'
        )
    return kwargs


def subplots(
    array=None, ncols=1, nrows=1,
    ref=1, order='C',
    aspect=1, figsize=None,
    width=None, height=None, journal=None,
    axwidth=None, axheight=None,
    hspace=None, wspace=None, space=None,
    hratios=None, wratios=None,
    width_ratios=None, height_ratios=None,
    left=None, bottom=None, right=None, top=None,
    basemap=None, proj=None, projection=None,
    proj_kw=None, projection_kw=None,
    **kwargs
):
    """
    Create a figure with a single subplot or arbitrary grids of subplots,
    analogous to `matplotlib.pyplot.subplots`. The subplots can be drawn with
    arbitrary projections.

    Parameters
    ----------
    array : 2d array-like of int, optional
        Array specifying complex grid of subplots. Think of
        this array as a "picture" of your figure. For example, the array
        ``[[1, 1], [2, 3]]`` creates one long subplot in the top row, two
        smaller subplots in the bottom row. Integers must range from 1 to the
        number of plots.

        ``0`` indicates an empty space. For example, ``[[1, 1, 1], [2, 0, 3]]``
        creates one long subplot in the top row with two subplots in the bottom
        row separated by a space.
    ncols, nrows : int, optional
        Number of columns, rows. Ignored if `array` was passed.
        Use these arguments for simpler subplot grids.
    order : {'C', 'F'}, optional
        Whether subplots are numbered in column-major (``'C'``) or row-major
        (``'F'``) order. Analogous to `numpy.array` ordering. This controls
        the order that subplots appear in the `SubplotsContainer` returned by
        this function, and the order of subplot a-b-c labels (see
        `~proplot.axes.Axes.format`).
    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    width, height : float or str, optional
        The figure width and height. If you specify just one, the aspect
        ratio `aspect` of the reference subplot `ref` will be preserved.
    ref : int, optional
        The reference subplot number. The `axwidth`, `axheight`, and `aspect`
        keyword args are applied to this subplot, and the aspect ratio is
        conserved for this subplot in the tight layout adjustment. If you
        did not specify `width_ratios` and `height_ratios`, the `axwidth`,
        `axheight`, and `aspect` settings will apply to *all* subplots --
        not just the `ref` subplot.
    axwidth, axheight : float or str, optional
        The width, height of the reference subplot. Units are interpreted by
        `~proplot.utils.units`. Default is :rc:`subplots.axwidth`. Ignored
        if `width`, `height`, or `figsize` was passed.
    aspect : float or length-2 list of floats, optional
        The reference subplot aspect ratio, in numeric form (width divided by
        height) or as a (width, height) tuple. Ignored if both `width` *and*
        `height` or both `axwidth` *and* `axheight` were passed.
    width_ratios, height_ratios : float or list thereof, optional
        Passed to `~proplot.gridspec.GridSpec`, denotes the width
        and height ratios for the subplot grid. Length of `width_ratios`
        must match the number of rows, and length of `height_ratios` must
        match the number of columns.
    wratios, hratios
        Aliases for `width_ratios`, `height_ratios`.
    wspace, hspace, space : float or str or list thereof, optional
        Passed to `~proplot.gridspec.GridSpec`, denotes the
        spacing between grid columns, rows, and both, respectively. If float
        or string, expanded into lists of length ``ncols - 1`` (for `wspace`)
        or length ``nrows - 1`` (for `hspace`).

        Units are interpreted by `~proplot.utils.units` for each element of
        the list. By default, these are determined by the "tight
        layout" algorithm.
    left, right, top, bottom : float or str, optional
        Passed to `~proplot.gridspec.GridSpec`, denotes the width of padding
        between the subplots and the figure edge. Units are interpreted by
        `~proplot.utils.units`. By default, these are determined by the
        "tight layout" algorithm.
    proj, projection : str, `cartopy.crs.Projection`, `~mpl_toolkits.basemap.Basemap`, \
list thereof, or dict thereof, optional
        The map projection specification(s). If ``'cartesian'`` (the default), a
        `~proplot.axes.CartesianAxes` is created. If ``'polar'``, a
        `~proplot.axes.PolarAxes` is created. Otherwise, the argument is
        interpreted by `~proplot.constructor.Proj`, and the result is used
        to make a `~proplot.axes.GeoAxes` (in this case the argument can be
        a `cartopy.crs.Projection` instance, a `~mpl_toolkits.basemap.Basemap`
        instance, or a projection name listed in :ref:`this table <proj_table>`).

        To use different projections for different subplots, you have
        two options:

        * Pass a *list* of projection specifications, one for each subplot.
          For example, ``plot.subplots(ncols=2, proj=('cartesian', 'robin'))``.
        * Pass a *dictionary* of projection specifications, where the
          keys are integers or tuples of integers that indicate the projection
          to use for the corresponding subplot(s). If a key is not provided, the
          default projection ``'cartesian'`` is used. For example,
          ``plot.subplots(ncols=4, proj={2: 'merc', (3, 4): 'stere'})`` creates
          a figure with a Cartesian axes for the first subplot, a Mercator
          projection for the second subplot, and a Stereographic projection
          for the third and fourth subplots.

    proj_kw, projection_kw : dict, list of dict, or dict of dicts, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or cartopy
        `~cartopy.crs.Projection` classes on instantiation. If dictionary of
        properties, applies globally. If list of dictionaries or dictionary of
        dictionaries, these apply to specific subplots, as with `proj`. For example,
        ``plot.subplots(ncols=2, proj_kw={1: {'lon_0': 0}, 2: {'lon_0': 180}})``
        centers the projection in the left subplot on the prime meridian and
        in the right subplot on the international dateline.
    basemap : bool, list of bool, or dict of bool, optional
        Passed to `~proplot.constructor.Proj`, determines whether projection
        string names like ``'pcarree'`` are used to create `~proplot.axes.BasemapAxes`
        or `~proplot.axes.CartopyAxes`. Default is ``False``. If boolean, applies
        to all subplots. If list or dict, applies to specific subplots, as with `proj`.
    journal : str, optional
        String name corresponding to an academic journal standard that is used
        to control the figure width and, if specified, the height. See the
        below table.

        .. _journal_table:

        ===========  ====================  ===============================================================================
        Key          Size description      Organization
        ===========  ====================  ===============================================================================
        ``'aaas1'``  1-column              `American Association for the Advancement of Science <aaas_>`_ (e.g. *Science*)
        ``'aaas2'``  2-column              ”
        ``'agu1'``   1-column              `American Geophysical Union <agu_>`_
        ``'agu2'``   2-column              ”
        ``'agu3'``   full height 1-column  ”
        ``'agu4'``   full height 2-column  ”
        ``'ams1'``   1-column              `American Meteorological Society <ams_>`_
        ``'ams2'``   small 2-column        ”
        ``'ams3'``   medium 2-column       ”
        ``'ams4'``   full 2-column         ”
        ``'nat1'``   1-column              `Nature Research <nat_>`_
        ``'nat2'``   2-column              ”
        ``'pnas1'``  1-column              `Proceedings of the National Academy of Sciences <pnas_>`_
        ``'pnas2'``  2-column              ”
        ``'pnas3'``  landscape page        ”
        ===========  ====================  ===============================================================================

        .. _aaas: https://www.sciencemag.org/authors/instructions-preparing-initial-manuscript
        .. _agu: https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Graphic-Requirements
        .. _ams: https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
        .. _nat: https://www.nature.com/nature/for-authors/formatting-guide
        .. _pnas: https://www.pnas.org/page/authors/format


    Other parameters
    ----------------
    **kwargs
        Passed to `~proplot.figure.Figure`.

    Returns
    -------
    f : `~proplot.figure.Figure`
        The figure instance.
    axs : `SubplotsContainer`
        A special list of axes instances. See `SubplotsContainer`.
    """  # noqa
    # Build array
    if order not in ('C', 'F'):  # better error message
        raise ValueError(
            f'Invalid order {order!r}. Choose from "C" (row-major, default) '
            f'and "F" (column-major).'
        )
    if array is None:
        array = np.arange(1, nrows * ncols + 1)[..., None]
        array = array.reshape((nrows, ncols), order=order)
    # Standardize array
    try:
        array = np.array(array, dtype=int)  # enforce array type
        if array.ndim == 1:  # interpret as single row or column
            array = array[None, :] if order == 'C' else array[:, None]
        elif array.ndim != 2:
            raise ValueError(
                f'Array must be 1-2 dimensional, but got {array.ndim} dims.'
            )
        array[array == None] = 0  # use zero for placeholder  # noqa
    except (TypeError, ValueError):
        raise ValueError(
            f'Invalid subplot array {array!r}. '
            'Must be 1d or 2d array of integers.'
        )
    # Get other props
    nums = np.unique(array[array != 0])
    naxs = len(nums)
    if {*nums.flat} != {*range(1, naxs + 1)}:
        raise ValueError(
            f'Invalid subplot array {array!r}. Numbers must span integers '
            '1 to naxs (i.e. cannot skip over numbers), with 0 representing '
            'empty spaces.'
        )
    if ref not in nums:
        raise ValueError(
            f'Invalid reference number {ref!r}. For array {array!r}, must be '
            f'one of {nums}.'
        )
    nrows, ncols = array.shape

    # Get some axes properties, where locations are sorted by axes id.
    # NOTE: These ranges are endpoint exclusive, like a slice object!
    # NOTE: 0 stands for empty
    axids = [np.where(array == i) for i in np.sort(np.unique(array)) if i > 0]
    xrange = np.array([[x.min(), x.max()] for _, x in axids])
    yrange = np.array([[y.min(), y.max()] for y, _ in axids])
    xref = xrange[ref - 1, :]  # range for reference axes
    yref = yrange[ref - 1, :]

    # Get basemap.Basemap or cartopy.crs.Projection instances for map
    proj = _not_none(projection=projection, proj=proj)
    proj = _axes_dict(naxs, proj, kw=False, default='cartesian')
    proj_kw = _not_none(projection_kw=projection_kw, proj_kw=proj_kw) or {}
    proj_kw = _axes_dict(naxs, proj_kw, kw=True)
    basemap = _axes_dict(naxs, basemap, kw=False, default=None)
    axes_kw = {num: {} for num in range(1, naxs + 1)}  # store add_subplot args
    for num, name in proj.items():
        # The default is CartesianAxes
        if name is None or name == 'cartesian':
            axes_kw[num]['projection'] = 'cartesian'

        # Builtin matplotlib polar axes, just use my overridden version
        elif name == 'polar':
            axes_kw[num]['projection'] = 'polar2'
            if num == ref:
                aspect = 1

        # Custom Basemap and Cartopy axes
        else:
            m = constructor.Proj(name, basemap=basemap[num], **proj_kw[num])
            package = m._proj_package
            if num == ref:
                if package == 'basemap':
                    aspect = (m.urcrnrx - m.llcrnrx) / (m.urcrnry - m.llcrnry)
                else:
                    aspect = (np.diff(m.x_limits) / np.diff(m.y_limits))[0]
            axes_kw[num].update({'projection': package, 'map_projection': m})

    # Figure and/or axes dimensions
    names, values = (), ()
    if journal:
        # if user passed width=<string > , will use that journal size
        figsize = _journals(journal)
        spec = f'journal={journal!r}'
        names = ('axwidth', 'axheight', 'width')
        values = (axwidth, axheight, width)
        width, height = figsize
    elif figsize:
        spec = f'figsize={figsize!r}'
        names = ('axwidth', 'axheight', 'width', 'height')
        values = (axwidth, axheight, width, height)
        width, height = figsize
    elif width is not None or height is not None:
        spec = []
        if width is not None:
            spec.append(f'width={width!r}')
        if height is not None:
            spec.append(f'height={height!r}')
        spec = ', '.join(spec)
        names = ('axwidth', 'axheight')
        values = (axwidth, axheight)
    # Raise warning
    for name, value in zip(names, values):
        if value is not None:
            warnings._warn_proplot(
                f'You specified both {spec} and {name}={value!r}. '
                f'Ignoring {name!r}.'
            )

    # Standardized dimensions
    width, height = units(width), units(height)
    axwidth, axheight = units(axwidth), units(axheight)

    # Standardized user input border spaces
    left, right = units(left), units(right)
    bottom, top = units(bottom), units(top)

    # Standardized user input spaces
    wspace = np.atleast_1d(units(_not_none(wspace, space)))
    hspace = np.atleast_1d(units(_not_none(hspace, space)))
    if len(wspace) == 1:
        wspace = np.repeat(wspace, (ncols - 1,))
    if len(wspace) != ncols - 1:
        raise ValueError(
            f'Require {ncols-1} width spacings for {ncols} columns, '
            'got {len(wspace)}.'
        )
    if len(hspace) == 1:
        hspace = np.repeat(hspace, (nrows - 1,))
    if len(hspace) != nrows - 1:
        raise ValueError(
            f'Require {nrows-1} height spacings for {nrows} rows, '
            'got {len(hspace)}.'
        )

    # Standardized user input ratios
    wratios = np.atleast_1d(_not_none(
        width_ratios=width_ratios, wratios=wratios, default=1,
    ))
    hratios = np.atleast_1d(_not_none(
        height_ratios=height_ratios, hratios=hratios, default=1,
    ))
    if len(wratios) == 1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios) == 1:
        hratios = np.repeat(hratios, (nrows,))
    if len(wratios) != ncols:
        raise ValueError(f'Got {ncols} columns, but {len(wratios)} wratios.')
    if len(hratios) != nrows:
        raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')

    # Fill subplots_orig_kw with user input values
    # NOTE: 'Ratios' are only fixed for panel axes, but we store entire array
    wspace, hspace = wspace.tolist(), hspace.tolist()
    wratios, hratios = wratios.tolist(), hratios.tolist()
    subplots_orig_kw = {
        'left': left, 'right': right, 'top': top, 'bottom': bottom,
        'wspace': wspace, 'hspace': hspace,
    }

    # Apply default settings
    share = kwargs.get('share', None)
    sharex = _not_none(kwargs.get('sharex', None), share, rc['subplots.share'])
    sharey = _not_none(kwargs.get('sharey', None), share, rc['subplots.share'])

    left = _not_none(left, pgridspec._default_space('left'))
    right = _not_none(right, pgridspec._default_space('right'))
    bottom = _not_none(bottom, pgridspec._default_space('bottom'))
    top = _not_none(top, pgridspec._default_space('top'))

    wspace, hspace = np.array(wspace), np.array(hspace)  # also copies!
    wspace[wspace == None] = pgridspec._default_space('wspace', sharex)  # noqa: E711, E501
    hspace[hspace == None] = pgridspec._default_space('hspace', sharey)  # noqa: E711, E501

    wratios, hratios = list(wratios), list(hratios)
    wspace, hspace = list(wspace), list(hspace)

    # Parse arguments, fix dimensions in light of desired aspect ratio
    figsize, gridspec_kw, subplots_kw = pgridspec._calc_geometry(
        nrows=nrows, ncols=ncols,
        aspect=aspect, xref=xref, yref=yref,
        left=left, right=right, bottom=bottom, top=top,
        width=width, height=height, axwidth=axwidth, axheight=axheight,
        wratios=wratios, hratios=hratios, wspace=wspace, hspace=hspace,
        wpanels=[''] * ncols, hpanels=[''] * nrows,
    )
    fig = plt.figure(
        FigureClass=pfigure.Figure, figsize=figsize, ref=ref,
        gridspec_kw=gridspec_kw, subplots_kw=subplots_kw,
        subplots_orig_kw=subplots_orig_kw,
        **kwargs
    )
    gridspec = fig._gridspec_main

    # Draw main subplots
    axs = naxs * [None]  # list of axes
    for idx in range(naxs):
        # Get figure gridspec ranges
        num = idx + 1
        x0, x1 = xrange[idx, 0], xrange[idx, 1]
        y0, y1 = yrange[idx, 0], yrange[idx, 1]
        # Draw subplot
        subplotspec = gridspec[y0:y1 + 1, x0:x1 + 1]
        with fig._context_authorize_add_subplot():
            axs[idx] = fig.add_subplot(
                subplotspec, number=num, main=True,
                **axes_kw[num]
            )

    # Shared axes setup
    # TODO: Figure out how to defer this to drawtime in #50
    # For some reason just adding _auto_share_setup() to draw() doesn't work
    for ax in axs:
        ax._auto_share_setup()

    # Return figure and axes
    n = ncols if order == 'C' else nrows
    return fig, SubplotsContainer(axs, n=n, order=order)


class SubplotsContainer(list):
    """
    List subclass and pseudo-2d array used as a container for the
    axes returned by `subplots`. See `~SubplotsContainer.__getattr__`
    and `~SubplotsContainer.__getitem__` for details.
    """
    def __init__(self, objs, n=1, order='C'):
        """
        Parameters
        ----------
        objs : list-like
            1d iterable of `~proplot.axes.Axes` instances.
        n : int, optional
            The length of the fastest-moving dimension, i.e. the number of
            columns when `order` is ``'C'``, and the number of rows when
            `order` is ``'F'``. Used to treat lists as pseudo-2d arrays.
        order : {'C', 'F'}, optional
            Whether 1d indexing returns results in row-major (C-style) or
            column-major (Fortran-style) order, respectively. Used to treat
            lists as pseudo-2d arrays.
        """
        if not all(isinstance(obj, paxes.Axes) for obj in objs):
            raise ValueError(
                f'Axes grid must be filled with Axes instances, got {objs!r}.'
            )
        super().__init__(objs)
        self._n = n
        self._order = order
        self._shape = (len(self) // n, n)[::(1 if order == 'C' else -1)]

    def __repr__(self):
        return 'SubplotsContainer([' + ', '.join(str(ax) for ax in self) + '])'

    def __setitem__(self, key, value):  # noqa: U100
        """
        Raise an error. This enforces pseudo immutability.
        """
        raise LookupError('SubplotsContainer is immutable.')

    def __getitem__(self, key):
        """
        If an integer is passed, the item is returned. If a slice is passed,
        a `SubplotsContainer` of the items is returned. You can also use 2D
        indexing, and the corresponding axes in the `SubplotsContainer` will
        be chosen.

        Example
        -------
        >>> import proplot as plot
        >>> fig, axs = plot.subplots(nrows=3, ncols=3, colorbars='b', bstack=2)
        >>> axs[0]  # the subplot in the top-right corner
        >>> axs[3]  # the first subplot in the second row
        >>> axs[1,2]  # the subplot in the second row, third from the left
        >>> axs[:,0]  # the subplots in the first column
        """
        # Allow 2d specification
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]
        # Do not expand single slice to list of integers or we get recursion!
        # len() operator uses __getitem__!
        if not isinstance(key, tuple):
            axlist = isinstance(key, slice)
            objs = list.__getitem__(self, key)
        elif len(key) == 2:
            axlist = any(isinstance(ikey, slice) for ikey in key)
            # Expand keys
            keys = []
            order = self._order
            for i, ikey in enumerate(key):
                if (i == 1 and order == 'C') or (i == 0 and order != 'C'):
                    n = self._n
                else:
                    n = len(self) // self._n
                if isinstance(ikey, slice):
                    start, stop, step = ikey.start, ikey.stop, ikey.step
                    if start is None:
                        start = 0
                    elif start < 0:
                        start = n + start
                    if stop is None:
                        stop = n
                    elif stop < 0:
                        stop = n + stop
                    if step is None:
                        step = 1
                    ikeys = [*range(start, stop, step)]
                else:
                    if ikey < 0:
                        ikey = n + ikey
                    ikeys = [ikey]
                keys.append(ikeys)

            # Get index pairs and get objects
            # Note that in double for loop, right loop varies fastest, so
            # e.g. axs[:,:] delvers (0,0), (0,1), ..., (0,N), (1,0), ...
            # Remember for order == 'F', SubplotsContainer was sent a list
            # unfurled in column-major order, so we replicate row-major
            # indexing syntax by reversing the order of the keys.
            objs = []
            if self._order == 'C':
                idxs = [
                    key0 * self._n + key1 for key0 in keys[0]
                    for key1 in keys[1]
                ]
            else:
                idxs = [
                    key1 * self._n + key0 for key1 in keys[1]
                    for key0 in keys[0]
                ]
            for idx in idxs:
                objs.append(list.__getitem__(self, idx))
            if not axlist:  # objs will always be length 1
                objs = objs[0]
        else:
            raise IndexError

        # Return
        if axlist:
            return SubplotsContainer(objs)
        else:
            return objs

    def __getattr__(self, attr):
        """
        If the attribute is *callable*, return a dummy function that loops
        through each identically named method, calls them in succession, and
        returns a tuple of the results. This lets us call arbitrary methods
        on several axes at once. If the `SubplotsContainer` has length ``1``,
        the single result is returned. If the attribute is *not callable*,
        returns a tuple of attributes for every object in the list.

        Example
        -------
        >>> import proplot as plot
        >>> fig, axs = plot.subplots(nrows=2, ncols=2)
        >>> axs.format(...)  # calls "format" on all subplots
        >>> panels = axs.panel_axes('right')  # returns SubplotsContainer of panels
        >>> panels.format(...)  # calls "format" on all panels
        """
        if not self:
            raise AttributeError(
                f'Invalid attribute {attr!r}, axes grid {self!r} is empty.'
            )
        objs = tuple(getattr(ax, attr) for ax in self)  # may raise error

        # Objects
        if not any(callable(_) for _ in objs):
            if len(self) == 1:
                return objs[0]
            else:
                return objs
        # Methods
        # NOTE: Must manually copy docstring because help() cannot inherit it
        elif all(callable(_) for _ in objs):
            @functools.wraps(objs[0])
            def _iterator(*args, **kwargs):
                result = []
                for func in objs:
                    result.append(func(*args, **kwargs))
                if len(self) == 1:
                    return result[0]
                elif all(res is None for res in result):
                    return None
                elif all(isinstance(res, paxes.Axes) for res in result):
                    return SubplotsContainer(result, n=self._n, order=self._order)
                else:
                    return tuple(result)
            _iterator.__doc__ = inspect.getdoc(objs[0])
            return _iterator

        # Mixed
        raise AttributeError(f'Found mixed types for attribute {attr!r}.')

    @property
    def shape(self):
        """
        The "shape" of the 2d subplot grid assumed when performing 2d
        indexing.  For :ref:`complex subplot grids <ug_intro>`, where
        subplots may span contiguous rows and columns, this "shape" may be
        incorrect. In such cases, 1d indexing should always be used.
        """
        return self._shape


# Deprecations
subplot_grid = warnings._rename_obj('subplot_grid', SubplotsContainer)
