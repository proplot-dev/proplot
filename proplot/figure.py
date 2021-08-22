#!/usr/bin/env python3
"""
The figure class used for all ProPlot figures.
"""
import functools
import inspect
import os
from collections.abc import MutableSequence
from numbers import Integral

import matplotlib.axes as maxes
import matplotlib.figure as mfigure
import matplotlib.gridspec as mgridspec
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import numpy as np

from . import axes as paxes
from . import constructor
from . import gridspec as pgridspec
from .config import _parse_format, _translate_loc, rc, rc_matplotlib
from .internals import ic  # noqa: F401
from .internals import (
    _empty_context,
    _not_none,
    _pop_params,
    _snippet_manager,
    _state_context,
    _version_mpl,
    docstring,
    warnings,
)
from .utils import units

__all__ = [
    'Figure',
    'SubplotGrid',
    'SubplotsContainer'  # deprecated
]


# Translation of input projections to native axes projections
# NOTE: We work outside of matplotlib's built-in registration system to avoid conflicts
AXES_PROJS = {
    'cartesian': 'proplot_cartesian',
    'cart': 'proplot_cartesian',
    'polar': 'proplot_polar',
    'three': 'proplot_three',
    '3d': 'proplot_three',
}
# Preset figure widths or sizes based on academic journal recommendations
# NOTE: Please feel free to add to this!
JOURNAL_SIZES = {
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


# Figure docstring
_figure_docstring = """
refnum : int, optional
    The reference subplot number. The `refwidth`, `refheight`, and `refaspect`
    keyword args are applied to this subplot, and the aspect ratio is conserved
    for this subplot in the `~Figure.auto_layout`. The default is the first
    subplot created in the figure.
refaspect : float or length-2 list of floats, optional
    The reference subplot aspect ratio. If scalar, this indicates the width
    divided by height. If 2-tuple, indicates the (width, height). Ignored if
    both `figwidth` *and* `figheight` or both `refwidth` *and* `refheight`
    were passed.
refwidth, refheight : float or str, optional
    The width, height of the reference subplot. Default is :rc:`subplots.refwidth`.
    %(units.in)s
    Ignored if `figwidth`, `figheight`, or `figsize` was passed.
ref, aspect, axwidth, axheight
    Aliases for `refnum`, `refaspect`, `refwidth`, `refheight`.
    *These may be deprecated in a future release.*
figwidth, figheight : float or str, optional
    The figure width and height.
    %(units.in)s
    If you specify just one, the aspect ratio `refaspect` of the reference
    subplot will be preserved.
width, height
    Aliases for `figwidth`, `figheight`.
figsize : length-2 tuple, optional
    Tuple specifying the figure ``(width, height)``.
sharex, sharey, share \
: {0, False, 1, 'labels', 'labs', 2, 'limits', 'lims', 3, True, 4, 'all'}, optional
    The axis sharing "level" for the *x* axis, *y* axis, or both axes.
    Default is :rc:`subplots.share`. Options are as follows:

    * ``0`` or ``False``: No axis sharing. This also sets the default `spanx`
      and `spany` values to ``False``.
    * ``1`` or ``'labels'`` or ``'labs'``: Only draw axis labels on the bottommost
      row or leftmost column of subplots. Tick labels still appear on every subplot.
    * ``2`` or ``'limits'`` or ``'lims'``: As above but force the axis limits, scales,
      and tick locations to be identical. Tick labels still appear on every subplot.
    * ``3`` or ``True``: As above but only show the tick labels on the bottommost
      row and leftmost column of subplots.
    * ``4`` or ``'all'``: As above but also share the axis limits, scales, and
      tick locations between subplots not in the same row or column.

spanx, spany, span : bool or {0, 1}, optional
    Whether to use "spanning" axis labels for the *x* axis, *y* axis, or both
    axes. Default is ``False`` if `sharex`, `sharey`, or `share` are ``0`` or
    ``False``, :rc:`subplots.span` otherwise. When ``True``, a single, centered
    axis label is used for all axes with bottom and left edges in the same
    row or column.  This can considerably redundancy in your figure.

    "Spanning" labels integrate with "shared" axes. For example,
    for a 3-row, 3-column figure, with ``sharey > 1`` and ``spany=1``,
    your figure will have 1 ylabel instead of 9.
alignx, aligny, align : bool or {0, 1}, optional
    Whether to `"align" axis labels \
<https://matplotlib.org/stable/gallery/subplots_axes_and_figures/align_labels_demo.html>`__
    for the *x* axis, *y* axis, or both axes. Aligned labels always appear in the
    same row or column. This Only has an effect when `spanx`, `spany`, or `span`
    are ``False``. Default is :rc:`subplots.align`.
%(gridspec.shared)s
%(gridspec.scalar)s
tight : bool, optional
    Whether to have `~Figure.auto_layout` include tight layout adjustments.
    Default is :rc:`subplots.tight`. If you manually specified a spacing in
    the call to `~proplot.ui.subplots`, it will be used to override the tight
    layout spacing. For example, with ``left=1``, the left margin is set to
    1 em-width, while the remaining margin widths are calculated automatically.
%(gridspec.tight)s
includepanels : bool, optional
    Whether to include panels when aligning figure "super titles" along the top
    of the subplot grid and when aligning the `spanx` *x* axis labels and `spany`
    *y* axis labels along the sides of the subplot grid. Default is ``False``.
mathtext_fallback : bool or str, optional
    Apply this :rc:`mathtext.fallback` value when drawing the figure. If
    ``True`` or string, unavailable glyphs are replaced with a glyph from a
    fallback font (Computer Modern by default). Otherwise, they are replaced
    with the "¤" dummy character. For details see this `mathtext tutorial \
<https://matplotlib.org/stable/tutorials/text/mathtext.html#custom-fonts>`__.
journal : str, optional
    String corresponding to an academic journal standard used to control the figure
    width `figwidth` and, if specified, the figure height `figheight`. See the below
    table. Feel free to add to this table by submitting a pull request.

    .. _journal_table:

    ===========  ====================  \
===============================================================================
    Key          Size description      Organization
    ===========  ====================  \
===============================================================================
    ``'aaas1'``  1-column              \
`American Association for the Advancement of Science <aaas_>`_ (e.g. *Science*)
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
    ``'pnas1'``  1-column              \
`Proceedings of the National Academy of Sciences <pnas_>`_
    ``'pnas2'``  2-column              ”
    ``'pnas3'``  landscape page        ”
    ===========  ====================  \
===============================================================================

    .. _aaas: \
https://www.sciencemag.org/authors/instructions-preparing-initial-manuscript
    .. _agu: \
https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Graphic-Requirements
    .. _ams: \
https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
    .. _nat: \
https://www.nature.com/nature/for-authors/formatting-guide
    .. _pnas: \
https://www.pnas.org/page/authors/format
"""
_snippet_manager['figure.figure'] = _figure_docstring


# Multiple subplots
_subplots_params_docstring = """
array : 2d array-like of int, optional
    Array specifying complex grid of subplots. Think of this as a "picture"
    of your figure. For example, the array ``[[1, 1], [2, 3]]`` creates one
    long subplot in the top row, two smaller subplots in the bottom row.
    Integers must range from 1 to the number of plots, and ``0`` indicates an
    empty space. For example, ``[[1, 1, 1], [2, 0, 3]]`` creates one long subplot
    in the top row with two subplots in the bottom row separated by a space.
ncols, nrows : int, optional
    The number of columns, rows in the subplot grid. Ignored if `array` was
    passed. Use these arguments for simpler subplot grids.
order : {'C', 'F'}, optional
    Whether subplots are numbered in column-major (``'C'``) or row-major (``'F'``)
    order. Analogous to `numpy.array` ordering. This controls the order that
    subplots appear in the `SubplotGrid` returned by this function, and the order
    of subplot a-b-c labels (see `~proplot.axes.Axes.format`).
%(axes.proj)s

    To use different projections for different subplots, you have
    two options:

    * Pass a *list* of projection specifications, one for each subplot.
      For example, ``pplt.subplots(ncols=2, proj=('cart', 'robin'))``.
    * Pass a *dictionary* of projection specifications, where the
      keys are integers or tuples of integers that indicate the projection
      to use for the corresponding subplot number(s). If a key is not
      provided, the default projection ``'cartesian'`` is used. For example,
      ``pplt.subplots(ncols=4, proj={2: 'cyl', (3, 4): 'stere'})`` creates
      a figure with a default Cartesian axes for the first subplot, a Mercator
      projection for the second subplot, and a Stereographic projection
      for the third and fourth subplots.

%(axes.proj_kw)s
    If dictionary of properties, applies globally. If list or dictionary of
    dictionaries, applies to specific subplots, as with `proj`. For example,
    ``pplt.subplots(ncols=2, proj='cyl', proj_kw=({'lon_0': 0}, {'lon_0': 180})``
    centers the projection in the left subplot on the prime meridian and in the
    right subplot on the international dateline.
%(axes.basemap)s
    If boolean, applies to all subplots. If list or dict, applies to specific
    subplots, as with `proj`.
%(gridspec.shared)s
%(gridspec.vector)s
%(gridspec.tight)s
"""
_snippet_manager['figure.subplots_params'] = _subplots_params_docstring


# Composed subplots docstring
_subplots_docstring = """
Add an arbitrary grid of subplots to the figure.

Parameters
----------
%(figure.subplots_params)s

Other parameters
----------------
%(figure.figure)s
**kwargs
    Passed to `Figure.add_subplot`.

Returns
-------
axs : `SubplotGrid`
    The axes instances stored in a `SubplotGrid`.

See also
--------
proplot.ui.figure
proplot.ui.subplots
proplot.figure.Figure
matplotlib.figure.Figure
proplot.figure.SubplotGrid
proplot.axes.Axes
"""
_snippet_manager['figure.subplots'] = _subplots_docstring


# Single subplots
_subplot_docstring = """
Add a subplot axes to the figure.

Parameters
----------
*args : int, tuple, or `~matplotlib.gridspec.SubplotSpec`, optional
    The subplot location specifier. Your options are:

    * A single 3-digit integer argument specifying the number of rows,
      number of columns, and gridspec number (using row-major indexing).
    * Three positional arguments specifying the number of rows, number of
      columns, and gridspec number (int) or number range (2-tuple of int).
    * A `~matplotlib.gridspec.SubplotSpec` instance generated by indexing
      a ProPlot `~proplot.gridspec.GridSpec`.

    For integer input, the implied geometry must be compatible with the implied
    geometry from previous calls -- for example, ``fig.add_subplot(331)`` followed
    by ``fig.add_subplot(132)`` is valid because the 1 row of the second input can
    be tiled into the 3 rows of the the first input, but ``fig.add_subplot(232)``
    will raise an error because 2 rows cannot be tiled into 3 rows. For
    `~matplotlib.gridspec.SubplotSpec` input, the `~matplotlig.gridspec.SubplotSpec`
    must be derived from the `~proplot.gridspec.GridSpec` used in previous calls.

    These restrictions arise because we allocate a single,
    unique `~Figure.gridspec` for each figure.
number : int, optional
    The axes number used for a-b-c labeling. See `~proplot.axes.Axes.format`
    for details. By default this is incremented automatically based on the other
    subplots in the figure. Use ``0`` or ``False`` to ensure the subplot has no
    a-b-c label. Note the number corresponding to ``a`` is ``1``, not ``0``.
autoshare : bool, optional
    Whether to automatically share the *x* and *y* axes with subplots spanning the
    same rows and columns based on the figure-wide `sharex` and `sharey` settings.
    Default is ``True``. This has no effect if :rcraw:`subplots.share` is ``False``
    or if ``sharex=False`` or ``sharey=False`` were passed to the figure.
%(axes.proj)s
%(axes.proj_kw)s
%(axes.basemap)s

Other parameters
----------------
**kwargs
    Passed to `matplotlib.axes.Axes`.
"""
_snippet_manager['figure.subplot'] = _subplot_docstring


# Single axes
_axes_docstring = """
Add a non-subplot axes to the figure.

Parameters
----------
rect : 4-tuple of float
    The (left, bottom, width, height) dimensions of the axes in
    figure-relative coordinates.
%(axes.proj)s
%(axes.proj_kw)s
%(axes.basemap)s
"""
_snippet_manager['figure.axes'] = _axes_docstring


# Colorbar or legend panel docstring
_space_docstring = """
loc : str, optional
    The {name} location. Valid location keys are as follows.

    ===========  =====================
    Location     Valid keys
    ===========  =====================
    left edge    ``'l'``, ``'left'``
    right edge   ``'r'``, ``'right'``
    bottom edge  ``'b'``, ``'bottom'``
    top edge     ``'t'``, ``'top'``
    ===========  =====================

row, rows
    Aliases for `span` for {name}s on the left or right side.
col, cols
    Aliases for `span` for {name}s on the top or bottom side.
span : int or 2-tuple of int, optional
    Describes how the {name} spans rows and columns of subplots.
    For example, ``fig.{name}(loc='b', col=1)`` draws a {name}
    beneath the leftmost column of subplots, and
    ``fig.{name}(loc='b', cols=(1,2))`` draws a {name} beneath the
    left two columns of subplots. By default, the {name} will span
    all rows and columns.
space : float or str, optional
    The fixed space between the {name} and the subplot grid. Units are
    interpreted by `~proplot.utils.units`. When the tight layout algorithm
    is active for the figure, this is adjusted automatically using `pad`.
    Otherwise, a suitable default is selected.
pad : float or str, optional
    The tight layout padding between the subplot grid and the {name}.
    Default is :rc:`subplots.innerpad` for the first {name} and
    :rc:`subplots.panelpad` for subsequently stacked {name}s.
"""
_snippet_manager['figure.colorbar_space'] = _space_docstring.format(name='colorbar')
_snippet_manager['figure.legend_space'] = _space_docstring.format(name='legend')


# Save docstring
_save_docstring = """
Save the figure.

Parameters
----------
path : path-like, optional
    The file path. User paths are expanded with `os.path.expanduser`.
**kwargs
    Passed to `~matplotlib.figure.Figure.savefig`

See also
--------
Figure.save
Figure.savefig
matplotlib.figure.Figure.savefig
"""
_snippet_manager['figure.save'] = _save_docstring


def _get_journal_size(preset):
    """
    Return the width and height corresponding to the given preset.
    """
    value = JOURNAL_SIZES.get(preset, None)
    if value is None:
        raise ValueError(
            f'Unknown preset figure size specifier {preset!r}. '
            'Current options are: '
            + ', '.join(map(repr, JOURNAL_SIZES.keys()))
        )
    figwidth = figheight = None
    try:
        figwidth, figheight = value
    except (TypeError, ValueError):
        figwidth = value
    return figwidth, figheight


def _add_canvas_preprocessor(canvas, method):
    """
    Return a pre-processer that can be used to override instance-level
    canvas draw() and print_figure() methods. This applies tight layout
    and aspect ratio-conserving adjustments and aligns labels. Required
    so canvas methods instantiate renderers with the correct dimensions.
    """
    # NOTE: Renderer must be (1) initialized with the correct figure size or
    # (2) changed inplace during draw, but vector graphic renderers *cannot*
    # be changed inplace. So options include (1) monkey patch
    # canvas.get_width_height, overriding figure.get_size_inches, and exploit
    # the FigureCanvasAgg.get_renderer() implementation (because FigureCanvasAgg
    # queries the bbox directly rather than using get_width_height() so requires
    # workaround), (2) override bbox and bbox_inches as *properties* (but these
    # are really complicated, dangerous, and result in unnecessary extra draws),
    # or (3) simply override canvas draw methods. Our choice is #3.
    def _canvas_preprocess(self, *args, **kwargs):
        fig = self.figure  # update even if not stale! needed after saves
        func = getattr(type(self), method)  # the original method

        # Bail out if we are already adjusting layout
        # NOTE: The _is_adjusting check necessary when inserting new
        # gridspec rows or columns with the qt backend.
        # NOTE: Return value for macosx _draw is the renderer, for qt draw is
        # nothing, and for print_figure is some figure object, but this block
        # has never been invoked when calling print_figure.
        if fig._is_adjusting:
            if method == '_draw':  # macosx backend
                return fig._get_renderer()
            else:
                return

        # Adjust layout
        # NOTE: The authorized_context is needed because some backends disable
        # constrained layout or tight layout before printing the figure.
        # NOTE: *Critical* to not add print_figure renderer to the cache when the print
        # method (print_pdf, print_png, etc.) calls Figure.draw(). Otherwise have issues
        # where (1) figure size and/or figure bounds are incorrect after saving figure
        # *then* displaying it in qt or inline notebook backends, and (2) figure fails
        # to update correctly after successively modifying and displaying within inline
        # notebook backend (previously worked around this by forcing additional draw()
        # call in this function before proceeding with print_figure).
        ctx1 = fig._context_adjusting(cache=(method != 'print_figure'))
        ctx2 = fig._context_authorized()  # backends might call set_constrained_layout()
        ctx3 = rc.context(fig._mathtext_context)  # draw with figure-specific setting
        with ctx1, ctx2, ctx3:
            fig.auto_layout()
            return func(self, *args, **kwargs)

    # Add preprocessor
    setattr(canvas, method, _canvas_preprocess.__get__(canvas))
    return canvas


class Figure(mfigure.Figure):
    """
    The `~matplotlib.figure.Figure` subclass used by proplot.
    """
    # Shared error and warning messages
    _share_message = (
        'Axis sharing level can be 0 or False (share nothing), '
        "1 or 'labels' or 'labs' (share axis labels), "
        "2 or 'limits' or 'lims' (share axis limits and axis labels), "
        '3 or True (share axis limits, axis labels, and tick labels), '
        "or 4 or 'all' (share axis labels and tick labels in the same gridspec "
        'rows and columns and share axis limits across all subplots).'
    )
    _space_message = (
        'To set the left, right, bottom, top, wspace, or hspace gridspec values, '
        'pass them as keyword arguments to pplt.subplots(). Please note they are '
        'now specified in physical units, with strings interpreted by pplt.units() '
        'and floats interpreted as font size-widths.'
    )
    _tight_message = (
        'ProPlot uses its own tight layout algorithm that is activated by default. '
        "To disable it, set pplt.rc['subplots.tight'] to False or pass tight=False "
        'to pplt.subplots(). For details, see fig.auto_layout().'
    )
    _warn_interactive = True  # disabled after first warning

    def __repr__(self):
        opts = {}
        for attr in ('refaspect', 'refwidth', 'refheight', 'figwidth', 'figheight'):
            value = getattr(self, '_' + attr)
            if value is not None:
                opts[attr] = np.round(value, 2)
        geom = ''
        if self.gridspec:
            nrows, ncols = self.gridspec.get_subplot_geometry()
            geom = f'nrows={nrows}, ncols={ncols}, '
        opts = ', '.join(f'{key}={value!r}' for key, value in opts.items())
        return f'Figure({geom}{opts})'

    # NOTE: If _rename_kwargs argument is an invalid identifier, it is
    # simply used in the warning message.
    @_snippet_manager
    @warnings._rename_kwargs(
        '0.7', axpad='innerpad', autoformat='pplt.rc.autoformat = {}'
    )
    def __init__(
        self, *, refnum=None, ref=None, refaspect=None, aspect=None,
        refwidth=None, refheight=None, axwidth=None, axheight=None,
        figwidth=None, figheight=None, width=None, height=None, journal=None,
        sharex=None, sharey=None, share=None,  # used for default spaces
        spanx=None, spany=None, span=None,
        alignx=None, aligny=None, align=None,
        left=None, right=None, top=None, bottom=None,
        wspace=None, hspace=None, space=None, wpad=None, hpad=None, pad=None,
        outerpad=None, innerpad=None, panelpad=None, includepanels=None,
        tight=None, mathtext_fallback=None,
        **kwargs
    ):
        """
        Parameters
        ----------
        %(figure.figure)s

        Other parameters
        ----------------
        **kwargs
            Passed to `matplotlib.figure.Figure`.

        See also
        --------
        proplot.ui.figure
        proplot.ui.subplots
        proplot.figure.Figure.add_subplot
        proplot.figure.Figure.subplots
        matplotlib.figure.Figure
        """
        # Add figure sizing settings
        # NOTE: We cannot catpure user-input 'figsize' here because it gets
        # automatically filled by the figure manager. See ui.figure().
        # NOTE: The figure size is adjusted according to these arguments by the
        # canvas preprocessor. Although in special case where both 'figwidth' and
        # 'figheight' were passes we update 'figsize' to limit side effects.
        refnum = _not_none(refnum=refnum, ref=ref, default=1)  # never None
        refaspect = _not_none(refaspect=refaspect, aspect=aspect)
        refwidth = _not_none(refwidth=refwidth, axwidth=axwidth)
        refheight = _not_none(refheight=refheight, axheight=axheight)
        figwidth = _not_none(figwidth=figwidth, width=width)
        figheight = _not_none(figheight=figheight, height=height)
        messages = []
        if journal is not None:
            jwidth, jheight = _get_journal_size(journal)
            if jwidth is not None and figwidth is not None:
                messages.append(('journal', journal, 'figwidth', figwidth))
            if jheight is not None and figheight is not None:
                messages.append(('journal', journal, 'figheight', figheight))
            figwidth = _not_none(jwidth, figwidth)
            figheight = _not_none(jheight, figheight)
        if figwidth is not None and refwidth is not None:
            messages.append(('figwidth', figwidth, 'refwidth', refwidth))
            refwidth = None
        if figheight is not None and refheight is not None:
            messages.append(('figheight', figheight, 'refheight', refheight))
            refheight = None
        if figwidth is None and figheight is None and refwidth is None and refheight is None:  # noqa: E501
            refwidth = rc['subplots.refwidth']  # always inches
        if np.iterable(refaspect):
            refaspect = refaspect[0] / refaspect[1]
        for key1, val1, key2, val2 in messages:
            warnings._warn_proplot(
                f'Got conflicting figure size arguments {key1}={val1!r} and '
                f'{key2}={val2!r}. Ignoring {key2!r}.'
            )
        self._refnum = refnum
        self._refaspect = refaspect
        self._refaspect_default = 1  # updated for imshow and geographic plots
        self._refwidth = units(refwidth, 'in')
        self._refheight = units(refheight, 'in')
        self._figwidth = units(figwidth, 'in')
        self._figheight = units(figheight, 'in')

        # Add special consideration for interactive backends
        backend = _not_none(rc.backend, '')
        backend = backend.lower()
        interactive = 'nbagg' in backend or 'ipympl' in backend
        if not interactive:
            pass
        elif figwidth is None or figheight is None:
            figsize = rc['figure.figsize']  # modified by proplot
            self._figwidth = figwidth = _not_none(figwidth, figsize[0])
            self._figheight = figheight = _not_none(figheight, figsize[1])
            self._refwidth = self._refheight = None  # critical!
            if self._warn_interactive:
                Figure._warn_interactive = False  # set class attribute
                warnings._warn_proplot(
                    'Auto-sized ProPlot figures are not compatible with interactive '
                    "backends like '%matplotlib widget' and '%matplotlib notebook'. "
                    f'Reverting to the figure size ({figwidth}, {figheight}). To make '
                    'auto-sized figures, please consider using the non-interactive '
                    '(default) backend. This warning message is shown the first time '
                    'you create a figure without explicitly specifying the size.'
                )

        # Add space settings
        # NOTE: This is analogous to 'subplotpars' but we don't worry about
        # user mutability. Think it's perfectly fine to ask users to simply
        # pass these to pplt.figure() or pplt.subplots(). Also overriding
        # 'subplots_adjust' would be confusing since we switch to absolute
        # units and that function is heavily used outside of proplot.
        params = {
            'left': left,
            'right': right,
            'top': top,
            'bottom': bottom,
            'wspace': wspace,
            'hspace': hspace,
            'space': space,
            'wpad': wpad,
            'hpad': hpad,
            'pad': pad,
            'outerpad': outerpad,
            'innerpad': innerpad,
            'panelpad': panelpad,
        }
        self._gridspec_params = params  # used to initialize the gridspec

        # Add tight layout setting and ignore native settings
        pars = kwargs.pop('subplotpars', None)
        if pars is not None:
            warnings._warn_proplot(
                f'Ignoring subplotpars={pars!r}. ' + self._space_message
            )
        if kwargs.pop('tight_layout', None):
            warnings._warn_proplot(
                'Ignoring tight_layout=True. ' + self._tight_message
            )
        if kwargs.pop('constrained_layout', None):
            warnings._warn_proplot(
                'Ignoring constrained_layout=True. ' + self._tight_message
            )
        if rc_matplotlib.get('figure.autolayout', False):
            warnings._warn_proplot(
                "Setting rc['figure.autolayout'] to False. " + self._tight_message
            )
        if rc_matplotlib.get('figure.constrained_layout.use', False):
            warnings._warn_proplot(
                "Setting rc['figure.constrained_layout.use'] to False. " + self._tight_message  # noqa: E501
            )
        rc_matplotlib['figure.autolayout'] = False  # this is rcParams
        rc_matplotlib['figure.constrained_layout.use'] = False  # this is rcParams
        self._autospace = _not_none(tight, rc['subplots.tight'])
        self._includepanels = _not_none(includepanels, False)

        # Translate share settings
        translate = {'labels': 1, 'labs': 1, 'limits': 2, 'lims': 2, 'all': 4}
        sharex = _not_none(sharex, share, rc['subplots.share'])
        sharey = _not_none(sharey, share, rc['subplots.share'])
        sharex = 3 if sharex is True else translate.get(sharex, sharex)
        sharey = 3 if sharey is True else translate.get(sharey, sharey)
        if sharex not in range(5):
            raise ValueError(f'Invalid sharex={sharex!r}. ' + self._share_message)
        if sharey not in range(5):
            raise ValueError(f'Invalid sharey={sharey!r}. ' + self._share_message)
        self._sharex = int(sharex)
        self._sharey = int(sharey)

        # Translate span and align settings
        spanx = _not_none(spanx, span, False if not sharex else None, rc['subplots.span'])  # noqa: E501
        spany = _not_none(spany, span, False if not sharey else None, rc['subplots.span'])  # noqa: E501
        if spanx and (alignx or align):  # only warn when explicitly requested
            warnings._warn_proplot('"alignx" has no effect when spanx=True.')
        if spany and (aligny or align):
            warnings._warn_proplot('"aligny" has no effect when spany=True.')
        self._spanx = bool(spanx)
        self._spany = bool(spany)
        alignx = _not_none(alignx, align, rc['subplots.align'])
        aligny = _not_none(aligny, align, rc['subplots.align'])
        self._alignx = bool(alignx)
        self._aligny = bool(aligny)

        # Initialize the figure
        # NOTE: Super labels are stored inside {axes: text} dictionaries
        self._gridspec = None
        self._panel_dict = {'left': [], 'right': [], 'bottom': [], 'top': []}
        self._subplot_dict = {}  # subplots indexed by number
        self._subplot_counter = 0  # avoid add_subplot() returning an existing subplot
        self._is_adjusting = False
        self._is_authorized = False
        with self._context_authorized():
            super().__init__(**kwargs)

        # Super labels. We don't rely on private matplotlib _suptitle attribute and
        # _align_axis_labels supports arbitrary spanning labels for subplot groups.
        # NOTE: Don't use 'anchor' rotation mode otherwise switching to horizontal
        # left and right super labels causes overlap. Current method is fine.
        self._suptitle = self.text(0.5, 0.95, '', ha='center', va='bottom')
        self._supxlabel_dict = {}  # an axes: label mapping
        self._supylabel_dict = {}  # an axes: label mapping
        self._suplabel_dict = {'left': {}, 'right': {}, 'bottom': {}, 'top': {}}
        self._suptitle_pad = rc['suptitle.pad']
        d = self._suplabel_props = {}  # store the super label props
        d['left'] = {'va': 'center', 'ha': 'right'}
        d['right'] = {'va': 'center', 'ha': 'left'}
        d['bottom'] = {'va': 'top', 'ha': 'center'}
        d['top'] = {'va': 'bottom', 'ha': 'center'}
        d = self._suplabel_pad = {}  # store the super label padding
        d['left'] = rc['leftlabel.pad']
        d['right'] = rc['rightlabel.pad']
        d['bottom'] = rc['bottomlabel.pad']
        d['top'] = rc['toplabel.pad']

        # Text drawing behavior
        if mathtext_fallback is None:
            context = {}
        elif _version_mpl >= 3.4:
            context = {'mathtext.fallback': mathtext_fallback if isinstance(mathtext_fallback, str) else 'cm' if mathtext_fallback else None}  # noqa: E501
        else:
            context = {'mathtext.fallback_to_cm': bool(mathtext_fallback)}
        self._mathtext_context = context
        self.format(rc_mode=1)  # initiate super label and title settings

    def _context_adjusting(self, cache=True):
        """
        Prevent re-running auto layout steps due to draws triggered by figure
        resizes. Otherwise can get infinite loops.
        """
        kw = {'_is_adjusting': True}
        if not cache:
            kw['_cachedRenderer'] = None  # temporarily ignore it
        return _state_context(self, **kw)

    def _context_authorized(self):
        """
        Prevent warning message when internally calling no-op methods. Otherwise
        emit warnings to help new users.
        """
        return _state_context(self, _is_authorized=True)

    def _parse_proj(
        self, proj=None, projection=None, proj_kw=None, projection_kw=None,
        basemap=None, use_aspect=False, **kwargs
    ):
        """
        Translate the user-input projection into keyword arguments that can be passed
        to `~proplot.figure.Figure.add_subplot`. Also return a default aspect ratio.
        """
        # Parse arguments
        proj = _not_none(proj=proj, projection=projection, default='cartesian')
        proj_kw = _not_none(proj_kw=proj_kw, projection_kw=projection_kw, default={})
        if isinstance(proj, str):
            proj = proj.lower()
            proj = AXES_PROJS.get(proj, proj)

        # Redirect to a basemap or cartopy projection
        # NOTE: The default aspect should already be '1'.
        if proj in AXES_PROJS.values():
            name = proj
        else:
            m = constructor.Proj(
                proj, basemap=basemap, include_axes_projections=True, **proj_kw
            )
            if m._proj_package == 'basemap':
                aspect = (m.urcrnrx - m.llcrnrx) / (m.urcrnry - m.llcrnry)
            else:
                aspect, = np.diff(m.x_limits) / np.diff(m.y_limits)  # pull out of array
            name = 'proplot_' + m._proj_package  # attribute added by Proj()
            kwargs['map_projection'] = m
            if use_aspect:
                self._refaspect_default = aspect

        kwargs['projection'] = name
        return kwargs

    def _get_align_axes(self, side):
        """
        Return the main axes along the edge of the figure.
        """
        x, y = 'xy' if side in ('left', 'right') else 'yx'
        axs = self._subplot_dict.values()
        if not axs:
            return []
        ranges = np.array([ax._range_subplotspec(x) for ax in axs])
        edge = ranges[:, 0].min() if side in ('left', 'top') else ranges[:, 1].max()
        idx = 0 if side in ('left', 'top') else 1
        axs = [ax for ax in axs if ax._range_subplotspec(x)[idx] == edge]
        axs = [ax for ax in sorted(axs, key=lambda ax: ax._range_subplotspec(y)[0])]
        axs = [ax for ax in axs if ax.get_visible()]
        return axs

    def _get_align_coord(self, side, axs):
        """
        Return the figure coordinate for centering spanning axis labels or super titles.
        """
        # Get position in figure relative coordinates
        if not all(isinstance(ax, paxes.Axes) for ax in axs):
            raise RuntimeError('Axes must be proplot axes.')
        if not all(isinstance(ax, maxes.SubplotBase) for ax in axs):
            raise RuntimeError('Axes must be subplots.')
        s = 'y' if side in ('left', 'right') else 'x'
        axs = [ax._panel_parent or ax for ax in axs]  # deflect to main axes
        if self._includepanels:  # include panel short axes
            axs = [_ for ax in axs for _ in ax._iter_axes(panels=True, children=False)]
        ranges = np.array([ax._range_subplotspec(s) for ax in axs])
        min_, max_ = ranges[:, 0].min(), ranges[:, 1].max()
        ax_lo = axs[np.where(ranges[:, 0] == min_)[0][0]]
        ax_hi = axs[np.where(ranges[:, 1] == max_)[0][0]]
        box_lo = ax_lo.get_subplotspec().get_position(self)
        box_hi = ax_hi.get_subplotspec().get_position(self)
        if s == 'x':
            pos = 0.5 * (box_lo.x0 + box_hi.x1)
        else:
            pos = 0.5 * (box_lo.y1 + box_hi.y0)  # 'lo' is actually on top of figure
        ax = axs[(np.argmin(ranges[:, 0]) + np.argmax(ranges[:, 1])) // 2]
        other = getattr(ax, f'_share{s}')
        if other and other._panel_parent:  # deflect to shared panel axes
            ax = other
        return pos, ax

    def _get_offset_coord(self, side, axs, renderer, *, pad=None, extra=None):
        """
        Return the figure coordinate for offsetting super labels and super titles.
        """
        s = 'x' if side in ('left', 'right') else 'y'
        cs = []
        objs = tuple(_ for ax in axs for _ in ax._iter_axes(panels=True, children=True, hidden=True))  # noqa: E501
        objs = objs + (extra or ())  # e.g. top super labels
        for obj in objs:
            bbox = obj.get_tightbbox(renderer)  # cannot use cached bbox
            attr = s + 'max' if side in ('top', 'right') else s + 'min'
            c = getattr(bbox, attr)
            c = (c, 0) if side in ('left', 'right') else (0, c)
            c = self.transFigure.inverted().transform(c)
            c = c[0] if side in ('left', 'right') else c[1]
            cs.append(c)
        width, height = self.get_size_inches()
        if pad is None:
            pad = self._suplabel_pad[side] / 72
            pad = pad / width if side in ('left', 'right') else pad / height
        return min(cs) - pad if side in ('left', 'bottom') else max(cs) + pad

    def _get_renderer(self):
        """
        Get a renderer at all costs. See matplotlib's tight_layout.py.
        """
        if self._cachedRenderer:
            renderer = self._cachedRenderer
        else:
            canvas = self.canvas
            if canvas and hasattr(canvas, 'get_renderer'):
                renderer = canvas.get_renderer()
            else:
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                canvas = FigureCanvasAgg(self)
                renderer = canvas.get_renderer()
        return renderer

    def _add_axes_panel(self, ax, side=None, **kwargs):
        """
        Add an axes panel.
        """
        # Interpret args
        # NOTE: Axis sharing not implemented for figure panels, 99% of the
        # time this is just used as construct for adding global colorbars and
        # legends, really not worth implementing axis sharing
        if not isinstance(ax, paxes.Axes) or not isinstance(ax, maxes.SubplotBase):
            raise RuntimeError('Cannot add axes panels to a non-subplot axes.')
        ax = ax._panel_parent or ax  # redirect to main axes
        side = _translate_loc(side, 'panel', default='right')

        # Add and setup the panel accounting for index changes
        # NOTE: Always put tick labels on the 'outside'
        gs = self.gridspec
        if not gs:
            raise RuntimeError('The gridspec must be active.')
        ss, share = gs._insert_panel(side, ax, **kwargs)
        pax = self.add_subplot(ss, autoshare=False, number=False)
        pax._panel_side = side
        pax._panel_share = share
        pax._panel_parent = ax
        ax._panel_dict[side].append(pax)
        ax._auto_share()
        axis = pax.yaxis if side in ('left', 'right') else pax.xaxis
        getattr(axis, 'tick_' + side)()  # set tick and tick label position
        axis.set_label_position(side)  # set label position
        return pax

    def _add_figure_panel(
        self, side=None, span=None, row=None, col=None, rows=None, cols=None, **kwargs
    ):
        """
        Add a figure panel.
        """
        # Interpret args and enforce sensible keyword args
        side = _translate_loc(side, 'panel', default='right')
        if side in ('left', 'right'):
            for key, value in (('col', col), ('cols', cols)):
                if value is not None:
                    raise ValueError(f'Invalid keyword {key!r} for {side!r} panel.')
            span = _not_none(span=span, row=row, rows=rows)
        else:
            for key, value in (('row', row), ('rows', rows)):
                if value is not None:
                    raise ValueError(f'Invalid keyword {key!r} for {side!r} panel.')
            span = _not_none(span=span, col=col, cols=cols)

        # Add and setup panel
        # NOTE: This relies on panel slot obfuscation built into gridspec
        gs = self.gridspec
        if not gs:
            raise RuntimeError('The gridspec must be active.')
        ss, _ = gs._insert_panel(side, span, filled=True, **kwargs)
        pax = self.add_subplot(ss, autoshare=False, number=False)
        plist = self._panel_dict[side]
        plist.append(pax)
        pax._panel_side = side
        pax._panel_share = False
        pax._panel_parent = None
        return pax

    def _align_axis_label(self, x):
        """
        Align *x* and *y* axis labels in the perpendicular and parallel directions.
        """
        # NOTE: Always use 'align' if 'span' is True to get correct offset
        # NOTE: Must trigger axis sharing here so that super label alignment
        # with tight=False is valid. Kind of kludgey but oh well.
        seen = set()
        span = getattr(self, '_span' + x)
        align = getattr(self, '_align' + x)
        for ax in self._subplot_dict.values():
            if isinstance(ax, paxes.CartesianAxes):
                ax._apply_axis_sharing()  # always!
            else:
                continue
            pos = getattr(ax, x + 'axis').get_label_position()
            if ax in seen or pos not in ('bottom', 'left'):
                continue  # already aligned or cannot align
            axs = ax._get_span_axes(pos, panels=False)  # returns panel or main axes
            if len(axs) == 1 or any(getattr(ax, '_share' + x) for ax in axs):
                continue  # nothing to align or axes have parents
            seen.update(axs)
            if span or align:
                if hasattr(self, '_align_label_groups'):
                    group = self._align_label_groups[x]
                else:
                    group = getattr(self, '_align_' + x + 'label_grp', None)
                if group is not None:  # fail silently to avoid fragile API changes
                    for ax in axs[1:]:
                        group.join(axs[0], ax)  # add to grouper
            if span:
                self._update_axis_label(pos, axs)

    def _align_super_labels(self, side, renderer):
        """
        Adjust the position of super labels.
        """
        # NOTE: Ensure title is offset only here.
        for ax in self._subplot_dict.values():
            ax._apply_title_above()
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid side {side!r}.')
        labels = self._suplabel_dict[side]
        axs = tuple(ax for ax, label in labels.items() if label.get_text())
        if not axs:
            return
        c = self._get_offset_coord(side, axs, renderer)
        for label in labels.values():
            s = 'x' if side in ('left', 'right') else 'y'
            label.update({s: c})

    def _align_super_title(self, renderer):
        """
        Adjust the position of the super title.
        """
        if not self._suptitle.get_text():
            return
        axs = self._get_align_axes('top')  # returns outermost panels
        if not axs:
            return
        labs = tuple(t for t in self._suplabel_dict['top'].values() if t.get_text())
        pad = (self._suptitle_pad / 72) / self.get_size_inches()[1]
        x = self._get_align_coord('top', axs)[0]
        y = self._get_offset_coord('top', axs, renderer, pad=pad, extra=labs)
        self._suptitle.set_ha('center')
        self._suptitle.set_va('bottom')
        self._suptitle.set_position((x, y))

    def _update_axis_label(self, side, axs):
        """
        Update the aligned axis label for the input axes.
        """
        # NOTE: Previously we secretly used matplotlib axis labels for spanning labels,
        # offsetting them between two subplots if necessary. Now we track designated
        # 'super' labels and replace the actual labels with spaces so they still impact
        # the tight bounding box and thus allocate space for the spanning label.
        x, y = 'xy' if side in ('bottom', 'top') else 'yx'
        labs = getattr(self, '_sup' + x + 'label_dict')  # dict of spanning labels
        setpos = getattr(mtext.Text, 'set_' + y)
        axislist = [getattr(ax, x + 'axis') for ax in axs]

        # Get the central label and "super" label for parallel alignment.
        # Initialize the super label if one does not already exist.
        c, ax = self._get_align_coord(side, axs)  # returns panel axes
        axis = getattr(ax, x + 'axis')  # use the central axis
        label = labs.get(ax, None)
        if label is None and not axis.label.get_text().strip():
            return  # nothing to transfer from the normal label
        if label is not None and not label.get_text().strip():
            return  # nothing to update on the super label
        if label is None:
            props = ('ha', 'va', 'rotation', 'rotation_mode')
            label = labs[ax] = self.text(0, 0, '')
            label.update({p: getattr(axis.label, 'get_' + p)() for p in props})

        # Copy text from central label to spanning label
        # NOTE: Must use spaces rather than newlines, otherwise tight layout
        # won't make room. Reason is Text implementation (see Text._get_layout())
        ax._transfer_text(axis.label, label)  # text, color, and font properties
        space = '\n'.join(' ' * (1 + label.get_text().count('\n')))
        for axis in axislist:  # should include original 'axis'
            axis.label.set_text(space)

        # Update spanning label position
        t = mtransforms.IdentityTransform()  # set in pixels
        cx, cy = axis.label.get_position()
        if x == 'x':
            trans = mtransforms.blended_transform_factory(self.transFigure, t)
            coord = (c, cy)
        else:
            trans = mtransforms.blended_transform_factory(t, self.transFigure)
            coord = (cx, c)
        label.set_transform(trans)
        label.set_position(coord)

        # Add simple monkey patch to ensure positions stay in sync
        # NOTE: Simply using axis._update_label_position() when this is called
        # is not sufficient. Fails with e.g. inline backend.
        def _set_coord(self, *args, **kwargs):
            setpos(self, *args, **kwargs)
            setpos(label, *args, **kwargs)
        setattr(axis.label, 'set_' + y, _set_coord.__get__(axis.label))

    def _update_super_labels(self, side, labels, **kwargs):
        """
        Assign the figure super labels and update settings.
        """
        # Update the label parameters
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid side {side!r}.')
        kw = rc.fill(
            {
                'color': side + 'label.color',
                'rotation': side + 'label.rotation',
                'size': side + 'label.size',
                'weight': side + 'label.weight',
                'family': 'font.family'
            },
            context=True,
        )
        kw.update(kwargs)  # used when updating *existing* labels
        props = self._suplabel_props[side]
        props.update(kw)  # used when creating *new* labels

        # Get the label axes
        # WARNING: In case users added labels then changed the subplot geometry we
        # have to remove labels whose axes don't match the current 'align' axes.
        axs = self._get_align_axes(side)
        if not labels or not axs:
            return  # occurs if called while adding axes
        if len(labels) != len(axs):
            raise ValueError(
                f'Got {len(labels)} {side} labels but found {len(axs)} axes '
                f'along the {side} side of the figure.'
            )
        src = self._suplabel_dict[side]
        extra = src.keys() - set(axs)
        for ax in extra:  # e.g. while adding axes
            text = src[ax].get_text()
            if text:
                warnings._warn_proplot(f'Removing {side} label with text {text!r} from axes {ax.number}.')  # noqa: E501
            src[ax].remove()  # remove from the figure

        # Update the label text
        for ax, label in zip(axs, labels):
            if ax in src:
                obj = src[ax]
            elif side in ('left', 'right'):
                trans = mtransforms.blended_transform_factory(self.transFigure, ax.transAxes)  # noqa: E501
                obj = src[ax] = self.text(0, 0.5, '', transform=trans)
                obj.update(props)
            else:
                trans = mtransforms.blended_transform_factory(ax.transAxes, self.transFigure)  # noqa: E501
                obj = src[ax] = self.text(0.5, 0, '', transform=trans)
                obj.update(props)
            if kw:
                obj.update(kw)
            if label is not None:
                obj.set_text(label)

    def _update_super_title(self, title, **kwargs):
        """
        Assign the figure super title and update settings.
        """
        kw = rc.fill(
            {
                'size': 'suptitle.size',
                'weight': 'suptitle.weight',
                'color': 'suptitle.color',
                'family': 'font.family'
            },
            context=True,
        )
        kw.update(kwargs)
        if kw:
            self._suptitle.update(kw)
        if title is not None:
            self._suptitle.set_text(title)

    @docstring._concatenate_original
    @_snippet_manager
    def add_axes(self, rect, **kwargs):
        """
        %(figure.axes)s
        """
        kwargs = self._parse_proj(**kwargs)
        return super().add_axes(rect, **kwargs)

    @docstring._concatenate_original
    @_snippet_manager
    def add_subplot(self, *args, number=None, **kwargs):
        """
        %(figure.subplot)s
        """
        # Parse arguments
        # Also set the 'refaspect_default' is this is the reference axes
        kwargs['use_aspect'] = number == self._refnum
        kwargs = self._parse_proj(**kwargs)
        args = args or (1, 1, 1)
        gs = self.gridspec

        # Integer
        if len(args) == 1 and isinstance(args[0], Integral):
            if not 111 <= args[0] <= 999:
                raise ValueError(f'Input {args[0]} must fall between 111 and 999.')
            args = tuple(map(int, str(args[0])))

        # Subplot spec
        if (
            len(args) == 1
            and isinstance(args[0], (maxes.SubplotBase, mgridspec.SubplotSpec))
        ):
            ss = args[0]
            if isinstance(ss, maxes.SubplotBase):
                ss = ss.get_subplotspec()
            if gs is None:
                gs = ss.get_topmost_subplotspec().get_gridspec()
            if not isinstance(gs, pgridspec.GridSpec):
                raise ValueError('Subplotspec must be derived from a proplot.GridSpec.')
            if ss.get_topmost_subplotspec().get_gridspec() is not gs:
                raise ValueError('Subplotspec must be derived from the active figure gridspec.')  # noqa: E501

        # Row and column spec
        # TODO: How to pass spacing parameters to gridspec? Consider overriding
        # subplots adjust? Or require using gridspec manually?
        elif (
            len(args) == 3
            and all(isinstance(arg, Integral) for arg in args[:2])
            and all(isinstance(arg, Integral) for arg in np.atleast_1d(args[2]))
        ):
            nrows, ncols, num = args
            i, j = np.resize(num, 2)
            if gs is None:
                gs = pgridspec.GridSpec(nrows, ncols)
            orows, ocols = gs.get_subplot_geometry()
            if orows % nrows:
                raise ValueError(f'Input rows {nrows} do not divide the figure gridspec rows {orows}.')  # noqa: E501
            if ocols % ncols:
                raise ValueError(f'Input columns {ncols} do not divide the figure gridspec columns {ocols}.')  # noqa: E501
            if any(_ < 1 or _ > nrows * ncols for _ in (i, j)):
                raise ValueError(f'Input subplot indices must fall between 1 and {nrows * ncols}. Got {i} and {j}.')  # noqa: E501
            rowfact, colfact = orows // nrows, ocols // ncols
            irow, icol = divmod(i - 1, ncols)  # convert to zero-based
            jrow, jcol = divmod(j - 1, ncols)
            irow, icol = irow * rowfact, icol * colfact
            jrow, jcol = (jrow + 1) * rowfact - 1, (jcol + 1) * colfact - 1
            ss = gs[irow:jrow + 1, icol:jcol + 1]

        # Otherwise
        else:
            raise ValueError(f'Invalid add_subplot positional arguments {args!r}.')

        # Add the subplot
        # NOTE: Pass subplotspec as keyword arg for mpl >= 3.4 workaround
        # NOTE: Must assign unique label to each subplot or else subsequent calls
        # to add_subplot() in mpl < 3.4 may return an already-drawn subplot in the
        # wrong location due to gridspec override. Is against OO package design.
        if number is None:
            number = 1 + max(self._subplot_dict, default=0)
        if number:  # must be added for a-b-c labels
            kwargs['number'] = number
        self._subplot_counter += 1  # unique label for each subplot
        kwargs.setdefault('label', f'subplot_{self._subplot_counter}')
        ax = super().add_subplot(ss, _subplot_spec=ss, **kwargs)
        if number:
            self._subplot_dict[number] = ax
        self.gridspec = gs  # trigger layout adjustment

        return ax

    @_snippet_manager
    def subplot(self, *args, **kwargs):  # shorthand
        """
        %(figure.subplot)s
        """
        return self.add_subplot(*args, **kwargs)

    @_snippet_manager
    def add_subplots(
        self, array=None, *, ncols=1, nrows=1, order='C',
        proj=None, projection=None, proj_kw=None, projection_kw=None, basemap=None,
        **kwargs
    ):
        """
        %(figure.subplots)s
        """
        # Clunky helper function
        # TODO: Consider deprecating and asking users to use add_subplot()
        def _axes_dict(naxs, input, kw=False, default=None):
            # First build up dictionary
            # 1. 'string' or {1: 'string1', (2, 3): 'string2'}
            if not kw:
                if np.iterable(input) and not isinstance(input, (str, dict)):
                    value = {num + 1: item for num, item in enumerate(input)}
                elif not isinstance(input, dict):
                    value = {range(1, naxs + 1): input}
            # 2. {'prop': value} or {1: {'prop': value1}, (2, 3): {'prop': value2}}
            else:
                nested = [isinstance(_, dict) for _ in input.values()]
                if not any(nested):  # any([]) == False
                    value = {range(1, naxs + 1): input.copy()}
                elif not all(nested):
                    raise ValueError(f'Invalid input {value!r}.')
            # Unfurl keys that contain multiple axes numbers
            output = {}
            for nums, item in value.items():
                nums = np.atleast_1d(nums)
                for num in nums.flat:
                    output[num] = item.copy() if kw else item
            # Fill with default values
            for num in range(1, naxs + 1):
                if num not in output:
                    output[num] = {} if kw else default
            if output.keys() != set(range(1, naxs + 1)):
                raise ValueError(
                    f'Have {naxs} axes, but {value!r} includes props for the axes: '
                    + ', '.join(map(repr, sorted(output))) + '.'
                )
            return output

        # Warning messages
        for key in ('gridspec_kw', 'subplot_kw'):
            kw = kwargs.pop(key, None)
            if not kw:
                continue
            warnings._warn_proplot(
                f'{key!r} is not necessary in ProPlot. Pass the '
                'parameters as keyword arguments instead.'
            )
            kwargs.update(kw or {})

        # Build subplot array
        if order not in ('C', 'F'):  # better error message
            raise ValueError(f"Invalid order={order!r}. Options are 'C' or 'F'.")
        if array is None:
            array = np.arange(1, nrows * ncols + 1)[..., None]
            array = array.reshape((nrows, ncols), order=order)
        array = np.atleast_1d(array)
        array[array == None] = 0  # None or 0 both valid placeholders  # noqa: E711
        array = array.astype(np.int)
        if array.ndim == 1:  # interpret as single row or column
            array = array[None, :] if order == 'C' else array[:, None]
        elif array.ndim != 2:
            raise ValueError(f'Expected 1d or 2d array of integers. Got {array}.')
        nums = np.unique(array[array != 0])
        naxs = len(nums)
        if any(num < 0 or not isinstance(num, Integral) for num in nums.flat):
            raise ValueError(f'Invalid subplot array {array!r}. Must contain positive integers.')  # noqa: E501
        nrows, ncols = array.shape

        # Get projection arguments used to initialize the axes
        proj = _not_none(projection=projection, proj=proj)
        proj = _axes_dict(naxs, proj, kw=False, default='cartesian')
        proj_kw = _not_none(projection_kw=projection_kw, proj_kw=proj_kw) or {}
        proj_kw = _axes_dict(naxs, proj_kw, kw=True)
        basemap = _axes_dict(naxs, basemap, kw=False)
        axes_kw = {
            num: {'proj': proj[num], 'proj_kw': proj_kw[num], 'basemap': basemap[num]}
            for num in proj
        }

        # Create gridspec and add subplots with subplotspecs
        # NOTE: The gridspec is added to the figure when we pass the subplotspec
        gridspec_kw = _pop_params(kwargs, pgridspec.GridSpec._update_params)
        gs = pgridspec.GridSpec(nrows, ncols, **gridspec_kw)
        axs = naxs * [None]  # list of axes
        axids = [np.where(array == i) for i in np.sort(np.unique(array)) if i > 0]
        axcols = np.array([[x.min(), x.max()] for _, x in axids])
        axrows = np.array([[y.min(), y.max()] for y, _ in axids])
        for idx in range(naxs):
            num = idx + 1
            x0, x1 = axcols[idx, 0], axcols[idx, 1]
            y0, y1 = axrows[idx, 0], axrows[idx, 1]
            ss = gs[y0:y1 + 1, x0:x1 + 1]
            kw = {**kwargs, **axes_kw[num], 'number': num}
            axs[idx] = self.add_subplot(ss, **kw)

        return SubplotGrid(axs)

    @_snippet_manager
    def subplots(self, *args, **kwargs):  # shorthand
        """
        %(figure.subplots)s
        """
        return self.add_subplots(*args, **kwargs)

    def auto_layout(self, renderer=None, aspect=None, tight=None, resize=None):
        """
        Automatically adjust the figure size and subplot positions. This is
        triggered automatically whenever the figure is drawn.

        Parameters
        ----------
        renderer : `~matplotlib.backend_bases.RendererBase`, optional
            The renderer. If ``None`` a default renderer will be produced.
        aspect : bool, optional
            Whether to update the figure size based on the reference subplot aspect
            ratio. By default, this is ``True``. This only has an effect if the
            aspect ratio is fixed (e.g., due to an image plot or geographic projection).
        tight : bool, optional
            Whether to update the figuer size and subplot positions according to
            a "tight layout". By default, this takes on the value of `tight` passed
            to `Figure`. If nothing was passed, it is :rc:`subplots.tight`.
        resize : bool, optional
            If ``False``, the current figure dimensions are fixed and automatic
            figure resizing is disabled. By default, the figure size may change
            unless both `figwidth` and `figheight` or `figsize` were passed
            to `~Figure.subplots`, `~Figure.set_size_inches` was called manually,
            or the figure was resized manually with an interactive backend.
        """
        # *Impossible* to get notebook backend to work with auto resizing so we
        # just do the tight layout adjustments and skip resizing.
        gs = self.gridspec
        renderer = self._get_renderer()
        if aspect is None:
            aspect = True
        if tight is None:
            tight = self._autospace
        if resize is False:  # fix the size
            self._figwidth, self._figheight = self.get_size_inches()
            self._refwidth = self._refheight = None  # critical!

        # Helper functions
        # NOTE: Have to draw legends and colorbars early (before reaching axes
        # draw methods) because we have to take them into account for alignment.
        # Also requires another figure resize (which triggers a gridspec update).
        def _draw_content():
            for ax in self._iter_axes(hidden=False, children=True):
                ax._draw_guides()  # may trigger resizes if panels are added
        def _align_content():  # noqa: E306
            for axis in 'xy':
                self._align_axis_label(axis)
            for side in ('left', 'right', 'top', 'bottom'):
                self._align_super_labels(side, renderer)
            self._align_super_title(renderer)

        # Update the layout
        # WARNING: Tried to avoid two figure resizes but made
        # subsequent tight layout really weird. Have to resize twice.
        _draw_content()
        if not gs:
            return
        if aspect:
            gs._auto_layout_aspect()
        _align_content()
        if tight:
            gs._auto_layout_space(renderer)
        _align_content()

    def format(
        self, axs=None,
        figtitle=None, suptitle=None, suptitle_kw=None,
        llabels=None, leftlabels=None, leftlabels_kw=None,
        rlabels=None, rightlabels=None, rightlabels_kw=None,
        blabels=None, bottomlabels=None, bottomlabels_kw=None,
        tlabels=None, toplabels=None, toplabels_kw=None,
        rowlabels=None, collabels=None, **kwargs,  # aliases
    ):
        """
        Call the ``format`` command for the input axes or
        for every subplot in the figure.

        Parameters
        ----------
        axs : list of axes, optional
            The list of axes. By default every subplot is used.
        %()

        Other parameters
        ----------------
        %(axes.rc)s
        **kwargs
            Passed to the ``format`` command of each axes.

        See also
        --------
        proplot.axes.Axes.format
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        """
        # Initiate context block
        rc_kw, rc_mode, kwargs = _parse_format(**kwargs)
        skip_axes = kwargs.pop('skip_axes', False)  # internal keyword arg
        with rc.context(rc_kw, mode=rc_mode):
            # Update background patch
            kw = rc.fill({'facecolor': 'figure.facecolor'}, context=True)
            self.patch.update(kw)

            # Update super title and label padding
            pad = rc.find('suptitle.pad', context=True)  # super title
            if pad is not None:
                self._suptitle_pad = pad
            for side in tuple(self._suplabel_pad):  # super labels
                pad = rc.find(side + 'label.pad', context=True)
                if pad is not None:
                    self._suplabel_pad[side] = pad

            # Update super title and labels text and settings
            suptitle_kw = suptitle_kw or {}
            leftlabels_kw = leftlabels_kw or {}
            rightlabels_kw = rightlabels_kw or {}
            bottomlabels_kw = bottomlabels_kw or {}
            toplabels_kw = toplabels_kw or {}
            self._update_super_title(
                _not_none(figtitle=figtitle, suptitle=suptitle),
                **suptitle_kw,
            )
            self._update_super_labels(
                'left',
                _not_none(rowlabels=rowlabels, leftlabels=leftlabels, llabels=llabels),
                **leftlabels_kw,
            )
            self._update_super_labels(
                'right',
                _not_none(rightlabels=rightlabels, rlabels=rlabels),
                **rightlabels_kw,
            )
            self._update_super_labels(
                'bottom',
                _not_none(bottomlabels=bottomlabels, blabels=blabels),
                **bottomlabels_kw,
            )
            self._update_super_labels(
                'top',
                _not_none(collabels=collabels, toplabels=toplabels, tlabels=tlabels),
                **toplabels_kw
            )

        # Update the main axes
        if skip_axes:  # avoid recursion
            return
        axs = axs or self._subplot_dict.values()
        for ax in axs:
            ax.format(rc_kw=rc_kw, rc_mode=rc_mode, skip_figure=True, **kwargs)

    @docstring._concatenate_original
    @_snippet_manager
    def colorbar(
        self, mappable, values=None, *, loc=None, location=None,
        row=None, col=None, rows=None, cols=None, span=None,
        space=None, pad=None, width=None, **kwargs
    ):
        """
        Draw a colorbar along the side of the figure.

        Parameters
        ----------
        %(axes.colorbar_args)s
        %(figure.colorbar_space)s
        length : float or str, optional
            The colorbar length. Units are relative to the span of the rows and
            columns of subplots. Default is :rc:`colorbar.length`.
        shrink : float, optional
            Alias for `length`. This is included for consistency with
            `matplotlib.figure.Figure.colorbar`.
        width : float or str, optional
            The colorbar width. Units are interpreted by
            `~proplot.utils.units`. Default is :rc:`colorbar.width`.

        Other parameters
        ----------------
        %(axes.colorbar_kwargs)s

        See also
        --------
        proplot.axes.Axes.colorbar
        matplotlib.figure.Figure.colorbar
        """
        ax = kwargs.pop('ax', None)
        cax = kwargs.pop('cax', None)
        # Fill this axes
        if cax is not None:
            with _state_context(cax, _internal_call=True):  # avoid wrapping pcolor
                return super().colorbar(mappable, cax=cax, **kwargs)
        # Axes panel colorbar
        elif ax is not None:
            return ax.colorbar(mappable, values, space=space, pad=pad, width=width, **kwargs)  # noqa: E501
        # Figure panel colorbar
        else:
            loc = _not_none(loc=loc, location=location, default='r')
            ax = self._add_figure_panel(loc, row=row, col=col, rows=rows, cols=cols, span=span, space=space, pad=pad, width=width)  # noqa: E501
            return ax.colorbar(mappable, values, loc='fill', **kwargs)

    @docstring._concatenate_original
    @_snippet_manager
    def legend(
        self, handles=None, labels=None, *, loc=None, location=None,
        row=None, col=None, rows=None, cols=None, span=None,
        space=None, pad=None, width=None, **kwargs
    ):
        """
        Draw a legend along the left, right, bottom, or top side of the
        figure, centered between the leftmost and rightmost (or topmost
        and bottommost) main subplots.

        Parameters
        ----------
        %(axes.legend_args)s
        %(figure.legend_space)s
        width : float or str, optional
            The space allocated for the legend box. This does nothing if the
            tight layout algorithm is active for the figure. Units are
            interpreted by `~proplot.utils.units`.

        Other parameters
        ----------------
        %(axes.legend_kwargs)s

        See also
        --------
        proplot.axes.Axes.legend
        matplotlib.axes.Axes.legend
        """
        ax = kwargs.pop('ax', None)
        # Axes panel legend
        if ax is not None:
            return ax.legend(handles, labels, space=space, pad=pad, width=width, **kwargs)  # noqa: E501
        # Figure panel legend
        else:
            loc = _not_none(loc=loc, location=location, default='r')
            ax = self._add_figure_panel(loc, row=row, col=col, rows=rows, cols=cols, span=span, space=space, pad=pad, width=width)  # noqa: E501
            return ax.legend(handles, labels, loc='fill', **kwargs)

    @_snippet_manager
    def save(self, filename, **kwargs):
        """
        %(figure.save)s
        """
        return self.savefig(filename, **kwargs)

    @docstring._concatenate_original
    @_snippet_manager
    def savefig(self, filename, **kwargs):
        """
        %(figure.save)s
        """
        # Automatically expand the user name. Undocumented because we
        # do not want to overwrite the matplotlib docstring.
        if isinstance(filename, str):
            filename = os.path.expanduser(filename)
        super().savefig(filename, **kwargs)

    @docstring._concatenate_original
    def set_canvas(self, canvas):
        """
        Set the figure canvas. Add monkey patches for the instance-level
        `~matplotlib.backend_bases.FigureCanvasBase.draw` and
        `~matplotlib.backend_bases.FigureCanvasBase.print_figure` methods.

        Parameters
        ----------
        canvas : `~matplotlib.backend_bases.FigureCanvasBase`
            The figure canvas.

        See also
        --------
        matplotlib.figure.Figure.set_canvas
        """
        # Set the canvas and add monkey patches to the instance-level draw and
        # print_figure methods. The latter is called by save() and by the inline
        # backend. See `_add_canvas_preprocessor` for details.
        _add_canvas_preprocessor(canvas, 'print_figure')
        if callable(getattr(canvas, '_draw', None)):  # for macosx backend
            _add_canvas_preprocessor(canvas, '_draw')
        else:
            _add_canvas_preprocessor(canvas, 'draw')
        super().set_canvas(canvas)

    def _is_same_size(self, figsize, eps=None):
        """
        Test if the figure size is unchanged up to some tolerance in inches.
        """
        eps = _not_none(eps, 0.01)
        figsize_active = self.get_size_inches()
        if figsize is None:  # e.g. GridSpec._calc_figsize() returned None
            return True
        else:
            return np.all(np.isclose(figsize, figsize_active, rtol=0, atol=eps))

    @docstring._concatenate_original
    def set_size_inches(self, w, h=None, *, forward=True, internal=False, eps=None):
        """
        Set the figure size. If this is being called manually or from an interactive
        backend, update the default layout with this fixed size. If the figure size is
        unchanged or this is an internal call, do not update the default layout.

        Parameters
        ----------
        *args : float
            The width and height passed as positional arguments or a 2-tuple.
        forward : bool, optional
            Whether to update the canvas.
        internal : bool, optional
            Whether this is an internal resize.
        eps : float, optional
            The deviation from the current size in inches required to treat this
            as a user-triggered figure resize that fixes the layout.

        See also
        --------
        matplotlib.figure.Figure.set_size_inches
        """
        # Parse input args
        figsize = w if h is None else (w, h)
        if not np.all(np.isfinite(figsize)):
            raise ValueError(f'Figure size must be finite, not {figsize}.')

        # Fix the figure size if this is a user action from an interactive backend
        # NOTE: If we fail to detect 'user' resize from the user, not only will
        # result be incorrect, but qt backend will crash because it detects a
        # recursive size change, since preprocessor size will differ.
        # NOTE: Bitmap renderers calculate the figure size in inches from
        # int(Figure.bbox.[width|height]) which rounds to whole pixels. When
        # renderer calls set_size_inches, size may be effectively the same, but
        # slightly changed due to roundoff error! Therefore only compare approx size.
        attrs = ('_is_idle_drawing', '_is_drawing', '_draw_pending')
        backend = any(getattr(self.canvas, attr, None) for attr in attrs)
        internal = internal or self._is_adjusting
        samesize = self._is_same_size(figsize, eps)
        context = _empty_context  # context not necessary most of the time
        if not backend and not internal and not samesize:
            context = self._context_adjusting  # do not trigger layout solver
            self._figwidth, self._figheight = figsize
            self._refwidth = self._refheight = None  # critical!

        # Apply the figure size
        # NOTE: If size changes we always update the gridspec to enforce fixed spaces
        # and panel widths (necessary since axes use figure relative coords)
        with context():  # avoid recursion
            super().set_size_inches(figsize, forward=forward)
        if not samesize:  # gridspec positions will resolve differently
            self.gridspec.update()

    def _iter_axes(self, hidden=False, children=False, panels=True):
        """
        Iterate over all axes and panels in the figure belonging to the
        `~proplot.axes.Axes` class. Exclude inset and twin axes.

        Parameters
        ----------
        hidden : bool, optional
            Whether to include "hidden" panels.
        children : bool, optional
            Whether to include child axes. Note this now includes "twin" axes.
        panels : bool or str or list of str, optional
            Whether to include panels or the panels to include.
        """
        # Parse panels
        if panels is False:
            panels = ()
        elif panels is True or panels is None:
            panels = ('left', 'right', 'bottom', 'top')
        elif isinstance(panels, str):
            panels = (panels,)
        if not set(panels) <= {'left', 'right', 'bottom', 'top'}:
            raise ValueError(f'Invalid sides {panels!r}.')
        # Iterate
        axs = (
            *self._subplot_dict.values(),
            *(ax for side in panels for ax in self._panel_dict[side])
        )
        for ax in axs:
            if not hidden and ax._panel_hidden:
                continue  # ignore hidden panel and its colorbar/legend child
            yield from ax._iter_axes(hidden=hidden, children=children, panels=panels)

    @property
    def gridspec(self):
        """
        The single `~proplot.gridspec.GridSpec` instance used for all
        subplots in the figure.
        """
        return self._gridspec

    @gridspec.setter
    def gridspec(self, gs):
        if not isinstance(gs, pgridspec.GridSpec):
            raise ValueError('Gridspec must be a proplot.GridSpec instance.')
        self._gridspec = gs
        gs.figure = self  # trigger copying settings from the figure

    @property
    def subplotgrid(self):
        """
        A `SubplotGrid` containing the numbered subplots in the figure,
        ordered by increasing subplot number.
        """
        return SubplotGrid([ax for num, ax in sorted(self._subplot_dict.items())])


# Add deprecated properties. There are *lots* of properties we pass to Figure
# and do not like idea of publicly tracking every single one of them. If we
# want to improve user introspection consider modifying Figure.__repr__.
for _attr in ('alignx', 'aligny', 'sharex', 'sharey', 'spanx', 'spany', 'tight', 'ref'):
    def _get_deprecated(self, attr=_attr):
        warnings._warn_proplot(
            f'The property {attr!r} is no longer public as of v0.8. It will be '
            'removed in a future release.'
        )
        return getattr(self, '_' + attr)
    _getter = property(_get_deprecated)
    setattr(Figure, _attr, property(_get_deprecated))


# Disable all native matplotlib layout and spacing functions
# Instead of using 'subplotpars' we record spacing params in LayoutSolver
# NOTE: In certain configurations with tight=False auto_layout() will not get
# called and thus LayoutSolver will not overwrite subplotpars. However it is not
# always clear when this is true. Better to disable user control of subplotpars.
for _attr, _msg in (
    ('set_tight_layout', Figure._tight_message),
    ('set_constrained_layout', Figure._tight_message),
    ('tight_layout', Figure._tight_message),
    ('init_layoutbox', Figure._tight_message),
    ('execute_constrained_layout', Figure._tight_message),
    ('subplots_adjust', Figure._space_message),
):
    _func = getattr(Figure, _attr, None)
    if _func is None:
        continue
    @functools.wraps(_func)  # noqa: E301
    def _disable_method(self, *args, func=_func, message=_msg, **kwargs):
        message = f'fig.{func.__name__}() has no effect on ProPlot figures. ' + message
        if self._is_authorized:
            return func(self, *args, **kwargs)
        else:
            warnings._warn_proplot(message)  # noqa: E501, U100
    _disable_method.__doc__ = None  # remove docs
    setattr(Figure, _attr, _disable_method)


class SubplotGrid(MutableSequence, list):
    """
    List-like object used to store subplots returned by `~Figure.subplots`. 1d indexing
    uses the underlying list of axes while 2d indexing uses the `~SubplotGrid.gridspec`.
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

    @docstring._obfuscate_signature  # hide deprecated args
    def __init__(self, iterable=None, n=None, order=None):
        """
        Parameters
        ----------
        iterable : list-like
            An iterable of `proplot.axes.Axes` subplots or their children.

        See also
        --------
        proplot.figure.Figure.subplots
        proplot.ui.subplots
        """
        if n is not None or order is not None:
            warnings._warn_proplot(
                f'Ignoring n={n!r} and order={order!r}. As of v0.8 SubplotGrid '
                'handles 2d indexing by leveraging the subplotspec extents rather than '
                'directly emulating 2d array indexing. These arguments are no longer '
                'needed and will be removed in a future release.'
            )
        iterable = _not_none(iterable, [])
        iterable = self._validate_item(iterable, scalar=False)
        super().__init__(iterable)

    def __getattr__(self, attr):
        """
        Get a missing attribute. Simply redirects to the axes if the `SubplotGrid`
        is singleton and raises an error otherwise. This can be very convenient
        for single-axes figures generated with `~Figure.subplots`.
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
        key : int, slice, or 2-tuple thereof
            The index. If 1d then the axes in the corresponding
            sublist are returned. If 2d then the axes that intersect
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
        elif isinstance(key, tuple) and len(key) == 2 and all(isinstance(ikey, (Integral, slice)) for ikey in key):  # noqa: E501
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
            The 1d index.
        value : `proplot.axes.Axes`
            The proplot subplot or its child or panel axes,
            or an iterable thereof if the index was a slice.
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
            'belonging to the same GridSpec.'
        )
        items = np.atleast_1d(items)
        if self:
            gridspec = self.gridspec  # compare against existing gridspec
        for item in items.flat:
            if not isinstance(item, paxes.Axes):
                raise ValueError(message)
            item = item._get_topmost_axes()
            if not isinstance(item, maxes.SubplotBase):
                raise ValueError(message)  # noqa: E501
            gs = item.get_subplotspec().get_gridspec()
            if not isinstance(gs, pgridspec.GridSpec) or (gridspec and gs is not gridspec):  # noqa: E501
                raise ValueError(message)
            gridspec = gs
        if not scalar:
            items = tuple(items.flat)
        elif items.size == 1:
            items = items.flat[0]
        else:
            raise ValueError('Input must be a single ProPlot axes.')
        return items

    @property
    def gridspec(self):
        """
        The `~proplot.gridspec.GridSpec` associated with the grid. This is used
        to resolve 2d indexing. See `~SubplotGrid.__getitem__` for details.
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
        # a 2d array-like object it should definitely have a shape attribute.
        return self.gridspec.get_subplot_geometry()


def _add_grid_command(name, command=None, seealso=None, returns_grid=False):
    # Build the docstring
    seealso = seealso or command or ()
    if isinstance(seealso, str):
        seealso = (seealso,)  # clean list of commands
    if command:
        command = f'`~{command}`'  # sphinx link to command
    else:
        command = f'``{name}``'  # literal text
    string = f"""
    Call {command} for every axes in the grid.

    Parameters
    ----------
    *args, **kwargs
        Passed to {command}.
    """
    string = inspect.cleandoc(string)
    string += '\n\nSee also\n--------\n' + '\n'.join(seealso)
    if returns_grid:
        string += '\n\nReturns\n-------\n'
        string += '    axs : `SubplotGrid`\n    A subplot grid of the results.'

    # Create the method
    def _grid_command(self, *args, **kwargs):
        objs = []
        for ax in self:
            obj = getattr(ax, name)(*args, **kwargs)
            if obj is not None:
                objs.append(obj)
        if not objs:
            return None
        elif all(isinstance(obj, paxes.Axes) for obj in objs):
            return SubplotGrid(objs)
        elif len(objs) == 1:
            return objs[0]
        else:
            return objs

    # Apply the method
    _grid_command.__name__ = name
    _grid_command.__doc__ = string
    setattr(SubplotGrid, name, _grid_command)


# Dynamically add commands
# TODO: Plot on axes in the grid along an extra input dimension?
_add_grid_command(
    'format',
    seealso=(f'proplot.axes.{s}Axes.format' for s in ('', 'Cartesian', 'Geo', 'Polar'))
)  # noqa: E501
for _name in (
    'panel',
    'panel_axes',
    'inset',
    'inset_axes',
    'altx',
    'alty',
    'dualx',
    'dualy',
    'twinx',
    'twiny',
    'text',
    'legend',
    'colorbar',
):
    if _name in ('altx', 'alty', 'twinx', 'twiny', 'dualx', 'dualy'):
        _command = f'proplot.axes.CartesianAxes.{_name}'
    else:
        _command = f'proplot.axes.Axes.{_name}'
    _returns_grid = _name not in ('text', 'legend', 'colorbar')
    _add_grid_command(_name, _command, returns_grid=_returns_grid)


# Deprecated
SubplotsContainer = warnings._rename_objs('0.8', SubplotsContainer=SubplotGrid)
