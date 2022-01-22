#!/usr/bin/env python3
"""
Utilities related to text labels.
"""
import matplotlib.offsetbox as moffsetbox
import matplotlib.patheffects as mpatheffects
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
from matplotlib import rcParams as rc_matplotlib

from . import ic  # noqa: F401
from . import _not_none, warnings

# Default border and box values
DEFAULT_BORDERCOLOR = 'w'
DEFAULT_BORDERINVERT = False
DEFAULT_BORDERWIDTH = 2
DEFAULT_BORDERSTYLE = 'miter'
DEFAULT_BOXFACECOLOR = 'w'
DEFAULT_BOXEDGECOLOR = 'k'
DEFAULT_BOXALPHA = 0.5
DEFAULT_BOXSTYLE = 'round'
DEFAULT_BOXPAD = 0.5


def _transfer_label(src, dest):
    """
    Transfer the input text object properties and content to the destination
    text object. Then clear the input object text.
    """
    if isinstance(src, moffsetbox.AnchoredText):
        src = src.txt
    elif not isinstance(src, mtext.Text):
        raise ValueError('Input must be Text or AnchoredText.')
    if isinstance(dest, moffsetbox.AnchoredText):
        dest = dest.txt._text
    elif not isinstance(dest, mtext.Text):
        raise ValueError('Input must be Text or AnchoredText.')
    text = src.get_text()
    dest.set_color(src.get_color())  # not a font property
    dest.set_fontproperties(src.get_fontproperties())  # size, weight, etc.
    if not text.strip():  # WARNING: must test strip() (see _align_axis_labels())
        return
    dest.set_text(text)
    src.set_text('')


@warnings._rename_kwargs(
    '0.10',
    bbox='box',
    bboxcolor='boxcolor',
    bboxalpha='boxalpha',
    bboxstyle='boxstyle',
    bboxpad='boxpad',
)
def _update_label(
    obj, border=None, box=None,
    bordercolor=None, borderwidth=None, borderinvert=None, borderstyle=None,
    boxcolor=None, boxalpha=None, boxstyle=None, boxpad=None, **kwargs
):
    """
    Update the text and (if applicable) offset box with "border"
    and "bbox" properties. This facillitates inset titles.
    """
    if isinstance(obj, mtext.Text):
        text = obj
        patch = obj.get_bbox_patch()  # NOTE: this can be None
    elif isinstance(obj, moffsetbox.AnchoredText):
        text = obj.txt._text
        patch = obj.patch
    else:
        raise ValueError('Input must be Text or AnchoredText.')
    text.update(kwargs)

    # Update text border
    border_props = {}
    if isinstance(border, dict):
        border_props.update(border)
        border = True
    if border is None:
        pass
    elif border:
        textcolor = text.get_color()
        patheffects = text.get_path_effects()
        if not patheffects:  # update with defaults
            bordercolor = _not_none(bordercolor, DEFAULT_BORDERCOLOR)
            borderinvert = _not_none(borderinvert, DEFAULT_BORDERINVERT)
            borderwidth = _not_none(borderwidth, DEFAULT_BORDERWIDTH)
            borderstyle = _not_none(borderstyle, DEFAULT_BORDERSTYLE)
        if borderinvert:
            bordercolor = _not_none(bordercolor, DEFAULT_BORDERCOLOR)
            textcolor, bordercolor = bordercolor, textcolor
        if borderwidth is not None:
            border_props.setdefault('linewidth', borderwidth)
        if bordercolor is not None:
            border_props.setdefault('foreground', bordercolor)
        if borderstyle is not None:
            border_props.setdefault('joinstyle', borderstyle)
        if textcolor is not None:
            text.set_color(textcolor)
        if patheffects:  # do not update with default
            patheffects[0].update(border_props)
        else:  # instantiate and apply defaults
            stroke = mpatheffects.Stroke(**border_props)
            text.set_path_effects([stroke, mpatheffects.Normal()])
    elif border is False:
        text.set_path_effects(None)

    # Update bounding box
    # NOTE: AnchoredOffsetbox padding is relative to rc['legend.fontsize']
    box_props = {}
    if isinstance(box, dict):
        box_props.update(box)
        box = True
    if box is None:
        pass
    elif box:
        edgecolor = None
        if patch is None or not patch.get_visible():
            edgecolor = DEFAULT_BOXEDGECOLOR
            boxcolor = _not_none(boxcolor, DEFAULT_BOXFACECOLOR)
            boxalpha = _not_none(boxalpha, DEFAULT_BOXALPHA)
            boxstyle = _not_none(boxstyle, DEFAULT_BOXSTYLE)
            boxpad = _not_none(boxpad, text.axes._title_pad)
        if edgecolor is not None:
            box_props.setdefault('edgecolor', edgecolor)
        if boxcolor is not None:
            box_props.setdefault('facecolor', boxcolor)
        if boxstyle is not None:
            box_props.setdefault('boxstyle', boxstyle)
        if boxalpha is not None:
            box_props.setdefault('alpha', boxalpha)
        if boxpad is not None:
            boxpad /= rc_matplotlib['legend.fontsize']  # convert points to em-widths
            box_props.setdefault('pad', boxpad)
        if patch is None:  # only possible for Text objects
            text.set_bbox(box_props)
        else:
            patch.set_visible(True)
            patch.update(box_props)
    else:
        if isinstance(obj, mtext.Text):
            text.set_bbox(None)  # disable the bbox
        else:
            patch.set_visible(False)
    return obj


class _AnchoredLabel(moffsetbox.AnchoredOffsetbox):
    """
    A class for storing anchored text. Embeds text in `~matplotlib.offsetbox.HPacker`
    instances for optional a-b-c label and title pair storage, supports adding and
    removing text from the packer, permits anchoring along bounding boxes of multiple
    axes (for `~proplot.gridspec.SubplotGrid.text` support), and allows automatically
    adjusted perpendicular offset and parallel alignment for "outer" locations when
    requested. The tight bounding box algorithm will include labels appearing
    between columns of the grid. A similar scheme might be used for auto-offsetting
    hanging twin axes and legends and removing the "panel" obfuscation.
    """
    def __init__(self, *args, loc=None, pad=None, sep=None, textprops=None, **kwargs):
        # Accept arbitrarily many texts
        from ..config import rc
        loc = _not_none(loc, rc['title.loc'])
        pad = _not_none(pad, rc['title.pad'])
        sep = _not_none(sep, rc['abc.titlepad'])
        txts = [moffsetbox.TextArea(arg, textprops=textprops) for arg in args]
        pack = moffsetbox.HPacker(pad=pad, sep=sep, children=txts)
        super().__init__(child=pack, **kwargs)
        self.patch.set_visible(False)

    def _update_offset_func(self, renderer, fontsize=None):
        # Update the offset function
        # TODO: Override this for only perpendicular offsets
        if fontsize is None:
            fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())

        def _offset(w, h, xd, yd, renderer):
            bbox = mtransforms.Bbox.from_bounds(0, 0, w, h)
            borderpad = self.borderpad * fontsize
            bbox_to_anchor = self.get_bbox_to_anchor()
            x0, y0 = self._get_anchored_bbox(self.loc, bbox, bbox_to_anchor, borderpad)
            self.axes._update_title_position(renderer)
            y1 = self.axes.bbox.y0
            _, y2 = self.axes.title.get_position()
            _, y2 = self.axes.title.get_transform().transform((0, y2))
            offset = y2 - y1
            return x0 + xd, y0 + yd + offset

        self.set_offset(_offset)

    def get_extent(self, renderer):
        # Get the inset extent after adjusting for the title
        # position and ignoring horizontal padding.
        # TODO: Override this for only perpendicular offsets
        w, h, xd, yd = self.get_child().get_extent(renderer)
        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        pad = self.pad * fontsize
        return w + 2 * pad, h + 2 * pad, xd + pad, yd + pad

    def get_coord(self, renderer):
        # Get the title offset coordinate axes
        # NOTE: Since proplot adds "twins" as child axes
        # they are already covered in the child axes iteration.
        title = self._child._children[1]
        x, _ = title.get_position()
        title.set_position((x, 1.0))
        axs = self._twinned_axes.get_siblings(self)
        for ax in self.child_axes:
            if ax is None:
                continue
            locator = ax.get_axes_locator()
            if locator:
                pos = locator(self, renderer)
                ax.apply_aspect(pos)
            else:
                ax.apply_aspect()
            axs = axs + [ax]
        # Get the coordinate
        # TODO: Align groups of labels by the same baseline. Possibly make label
        # groupings similar to twinned axes groupings. And possibly use custom
        # logic for "aligned" axis labels rather than using matplotlib logic.
        # NOTE: Since AnchoredOffsetbox already adjusts text baseline position by
        # its window extent we don't need to use the extra check employed in
        # matplotlib _update_title_position. Suggests this way is cleaner.
        top = 0
        for ax in axs:
            top = max(top, ax.bbox.ymax)
            if ax.xaxis.get_visible() and (
                ax.xaxis.get_label_position() == 'top'
                or ax.xaxis.get_ticks_position() in ('top', 'unknown')
            ):
                bb = ax.xaxis.get_tightbbox(renderer)
            else:
                bb = ax.get_window_extent(renderer)
            if bb is None:
                continue
            ymax = bb.ymax
            if ax.xaxis.get_visible() and any(
                tick.tick2line.get_visible() and not tick.label2.get_visible()
                for tick in ax.xaxis.majorTicks
            ):
                ymax += ax.xaxis.get_tick_padding()
            top = max(top, ymax)
        return top
