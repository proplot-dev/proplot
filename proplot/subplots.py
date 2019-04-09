#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures and axes.
The `subplots` function is all you'll need to directly use here.
It returns a `Figure` instance and an `axes_list` list of
`~proplot.axes.BaseAxes` axes.

Note that instead of separating various features into their own functions
(e.g.  a `generate_panel` function), we stuff most of them right into
the `subplots` function. The reason for this approach? We want to build a
**static "scaffolding"** before plotting anything, so that ProPlot can exert
a ton of control over the layout and make it look "nice" without any
manual tweaking on your part.

See `~Figure.smart_tight_layout` for details.
"""
import os
import re
import numpy as np
# Local modules, projection sand formatters and stuff
from .rcmod import rc
from .utils import _default, _timer, _counter, units, journals, ic
from . import axistools, gridspec, projs, axes
import functools
import warnings
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
import matplotlib.figure as mfigure
import matplotlib.transforms as mtransforms
# Aliases for panel names
_aliases = {
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'lpanel': 'leftpanel'
    }

#------------------------------------------------------------------------------#
# Miscellaneous stuff
#------------------------------------------------------------------------------#
# Wrapper functions, so user doesn't have to import pyplot
def close():
    """Alias for ``matplotlib.pyplot.close('all')``, included so you don't have
    to import `~matplotlib.pyplot`. Closes all figures stored
    in memory."""
    plt.close('all') # easy peasy

def show():
    """Alias for ``matplotlib.pyplot.show()``, included so you don't have
    to import `~matplotlib.pyplot`. Note this command should be
    unnecessary if you are doing inline iPython notebook plotting and ran the
    `~proplot.notebook.nbsetup` command."""
    plt.show()

# Helper classes
class axes_list(tuple):
    """Magical class that iterates through items and calls respective
    method (or retrieves respective attribute) on each one. For example,
    ``axs.format(color='r')`` colors all axes spines and tick marks red."""
    def __new__(cls, list_, order='C', n=1):
        """Adds static attribute to immutable object."""
        cls._order = order
        cls._n = n # means ncols or nrows, depending on order
        return super().__new__(cls, list_)

    def __repr__(self):
        """Wraps the string representation."""
        return 'axes_list(' + super().__repr__() + ')'

    def __getitem__(self, key):
        """Returns an `axes_list` version of the slice, or just the axes."""
        # Allow 2D specification
        # For weirder keys, raise error down the line
        # NOTE: When order=='F', panels were unfurled columnwise, so number
        # of rows is actually number of columns on the array
        # TODO: Test this stuff. For order=='F':
        # 0, 1 --> 1*nrows, but nrows is ncols
        # 2, 1 --> 1*nrows + 2
        if isinstance(key, tuple):
            if len(key)==1:
                key = key[0]
            elif len(key)==2 and all(isinstance(k, Number) for k in key):
                if self._order=='C':
                    key = key[0]*self._n + key[1]
                elif self._order=='F':
                    key = key[1]*self._n + key[0]
        axs = tuple.__getitem__(self, key)
        if isinstance(key, slice): # i.e. returns a list
            axs = axes_list(axs)
        return axs

    def __getattr__(self, attr):
        """Stealthily returns dummy function that actually loops through each
        axes method and calls them in succession. If attribute(s) are not
        callable, instead returns a tuple of the attributes."""
        # This is only called when __getattribute__ fails, so builtin
        # methods and hidden stuff like _n work fine
        attrs = *(getattr(ax, attr, None) for ax in self), # magical tuple expansion
        if None in attrs:
            raise AttributeError(f'Attribute "{attr}" not found.')
        elif all(callable(_) for _ in attrs):
            @functools.wraps(attrs[0])
            def iterator(*args, **kwargs):
                ret = []
                for attr in attrs:
                    ret.append(attr(*args, **kwargs))
                if len(ret)==1:
                    return ret[0]
                elif all(res is None for res in ret):
                    return None
                else:
                    return (*ret,) # expand to tuple
            return iterator
        elif not any(callable(_) for _ in attrs):
            # return attrs
            return axes_list(attrs)
        else:
            raise AttributeError(f'Found mixed types for attribute "{attr}".')

#------------------------------------------------------------------------------#
# Figure class
#------------------------------------------------------------------------------#
def _intervalx_errfix(ax, bbox):
    """Given an axes and a bounding box, pads the intervalx according to the
    matplotlib "tight layout" error associated with invisible y ticks."""
    if not isinstance(ax, axes.XYAxes):
        return bbox.intervalx
    xerr = ax._ytick_pad_error # error in x-direction, due to y ticks
    return (bbox.intervalx[0] + xerr[0], bbox.intervalx[1] - sum(xerr))

def _intervaly_errfix(ax, bbox):
    """Given an axes and a bounding box, pads the intervaly according to the
    matplotlib "tight layout" error associated with invisible x ticks."""
    if not isinstance(ax, axes.XYAxes):
        return bbox.intervaly
    yerr = ax._xtick_pad_error # error in y-direction, due to x ticks
    return (bbox.intervaly[0] + yerr[0], bbox.intervaly[1] - sum(yerr))

def _ax_span(ax, renderer, children=True):
    """Get span, accounting for panels, shared axes, and whether axes has
    been replaced by colorbar in same location."""
    # Get bounding boxes
    pairs = [(ax, ax.get_tightbbox(renderer))]
    ax = ax._colorbar_parent or ax
    if children:
        axs = (*ax.leftpanel, *ax.bottompanel, *ax.rightpanel, *ax.toppanel,
            ax._twinx_child, ax._twiny_child)
        for iax in axs:
            if not iax:
                continue
            iax = iax._colorbar_child or iax
            if iax.get_visible():
                pairs.append((iax, iax.get_tightbbox(renderer)))
    pairs = [(_,bbox) for _,bbox in pairs if bbox is not None]
    # Return arrays
    xs = np.array([_intervalx_errfix(ax,bbox) for ax,bbox in pairs])
    ys = np.array([_intervaly_errfix(ax,bbox) for ax,bbox in pairs])
    xspan = [xs[:,0].min(), xs[:,1].max()]
    yspan = [ys[:,0].min(), ys[:,1].max()]
    return xspan, yspan

def _ax_props(axs, renderer):
    """If this is a panel axes, check if user generated a colorbar axes,
    which stores all the artists/is what we really want to get a
    tight bounding box for."""
    if not axs: # e.g. an "empty panel" was passed
        return (np.empty((0,2)),)*4
    axs = [(ax._colorbar_child or ax) for ax in axs]
    axs = [ax for ax in axs if ax.get_visible()]
    if not axs:
        return (np.empty((0,2)),)*4
    spans = [_ax_span(ax, renderer) for ax in axs]
    xspans = np.array([span[0] for span in spans])
    yspans = np.array([span[1] for span in spans])
    xrange = np.array([(ax._colorbar_parent or ax)._xrange for ax in axs])
    yrange = np.array([(ax._colorbar_parent or ax)._yrange for ax in axs])
    return xspans, yspans, xrange, yrange

class Figure(mfigure.Figure):
    def __init__(self,
            tight=True, outertight=None, subplottight=None, paneltight=None,
            flush=False, wflush=None, hflush=None,
            outerpad=None, subplotpad=None, panelpad=None,
            rcreset=True, silent=True, # print various draw statements?
            **kwargs):
        """
        The `~matplotlib.figure.Figure` instance returned by `subplots`.

        Parameters
        ----------
        tight : None or bool, optional
            If not ``None``, overrides `outertight`, `subplottight`, and
            `paneltight`. Defaults to ``True``. These settings
            conveniently adjust subplot positioning and figure bounding box
            framing, without messing up subplot aspect ratios and panel widths.
            See `~Figure.smart_tight_layout` for details.
        outertight : None or bool, optional
            Whether to draw a tight bounding box around the whole figure.
            If ``None``, takes the value of `tight`.
        subplottight : None or bool, optional
            Whether to automatically space out subplots to prevent overlapping
            axis tick labels, etc. If ``None``, takes the value of `tight`.
        paneltight : None or bool, optional
            Whether to automatically space between subplots and their panels
            to prevent overlap. If ``None``, takes the value of `tight`.
        outerpad, subplotpad, panelpad : None, float, or str, optional
            Margin size for tight bounding box surrounding the edge of the
            figure, between subplots in the figure, and between panels and
            their parent subplots, respectively. If float, units are inches. If
            string, units are interpreted by `~proplot.utils.units`.
        flush, wflush, hflush : None or bool, optional
            Whether subplots should be "flush" against each other in the
            horizontal (`wflush`), vertical (`hflush`), or both (`flush`)
            directions. Useful if you want to let axis ticks overlap with other
            axes, and just want axes spines touching each other. Note that
            instead of ``flush=0``, you can also use ``subplottight=False``
            with manual gridspec spacings ``wspace=0`` and ``hspace=0``.
        rcreset : bool, optional
            Whether to reset all `~proplot.rcmod.rc` settings to their
            default values once the figure is drawn.
        silent : bool, optional
            Whether to silence messages related to drawing and printing
            the figure.
        **kwargs
            Passed to `matplotlib.figure.Figure`.
        """
        # NOTE: Do not document hidden args, or stuff that needs to
        # be set manually by subplots() after declaration.
        # Initialize figure with some custom attributes
        self._silent  = silent # whether to print message when we resize gridspec
        self._rcreset = rcreset
        # Tight toggling
        tight = _default(tight, rc['tight'])
        self._smart_tight_outer   = _default(outertight, tight)
        self._smart_tight_subplot = _default(subplottight, tight)
        self._smart_tight_panel   = _default(paneltight, tight)
        self._smart_tight_init = True # is figure in its initial state?
        # Padding and "flush" args
        self._extra_pad = 0 # sometimes matplotlib fails, cuts off super title! will add to this
        self._smart_outerpad = units(_default(outerpad, rc['subplot.outerpad']))
        self._smart_subplotpad  = units(_default(subplotpad,  rc['subplot.subplotpad']))
        self._smart_panelpad = units(_default(panelpad, rc['subplot.panelpad']))
        self._subplot_wflush = _default(flush, wflush)
        self._subplot_hflush = _default(flush, hflush)
        # Gridspec information
        # These are filled in by subplots()
        self._subplots_kw = None # extra special settings
        self._main_gridspec = None # gridspec encompassing drawing area
        self._main_axes = []  # list of 'main' axes (i.e. not insets or panels)
        self._span_labels = [] # add axis instances to this, and label position will be updated
        # Panels, initiate as empty
        self.leftpanel   = axes.EmptyPanel()
        self.bottompanel = axes.EmptyPanel()
        self.rightpanel  = axes.EmptyPanel()
        self.toppanel    = axes.EmptyPanel()
        # Initialize
        super().__init__(**kwargs)
        self.suptitle('') # this adds _suptitle position
        self._suptitle_transform = None

    def __getattribute__(self, attr, *args):
        """Enables the attribute aliases ``bpanel`` for ``bottompanel``,
        ``tpanel`` for ``toppanel``, ``lpanel`` for ``leftpanel``, and
        ``rpanel`` for ``rightpanel``."""
        attr = _aliases.get(attr, attr)
        return super().__getattribute__(attr, *args)

    def _lock_twins(self):
        """Lock shared axis limits to these axis limits. Used for secondary
        axis with alternate data scale."""
        # Helper func
        # Note negative heights should not break anything!
        def check(lim, scale):
            if re.match('^log', scale) and any(np.array(lim)<=0):
                raise ValueError('Axis limits go negative, and "alternate units" axis uses log scale.')
            elif re.match('^inverse', scale) and any(np.array(lim)<=0):
                raise ValueError('Axis limits cross zero, and "alternate units" axis uses inverse scale.')
        for ax in self._main_axes:
            # Match units for x-axis
            if ax._dualx_scale:
                # Get stuff
                # transform = mscale.scale_factory(twin.get_xscale(), twin.xaxis).get_transform()
                twin = ax._twiny_child
                offset, scale = ax._dualx_scale
                transform = twin.xaxis._scale.get_transform() # private API
                xlim_orig = ax.get_xlim()
                # Check, and set
                check(xlim_orig, twin.get_xscale())
                xlim = transform.inverted().transform(np.array(xlim_orig))
                if np.sign(np.diff(xlim_orig)) != np.sign(np.diff(xlim)): # the transform flipped it, so when we try to set limits, will get flipped again!
                    xlim = xlim[::-1]
                twin.set_xlim(offset + scale*xlim)
            # Match units for y-axis
            if ax._dualy_scale:
                # Get stuff
                # transform = mscale.scale_factory(twin.get_yscale(), ax.yaxis).get_transform()
                twin = ax._twinx_child
                offset, scale = ax._dualy_scale
                transform = twin.yaxis._scale.get_transform() # private API
                ylim_orig = ax.get_ylim()
                check(ylim_orig, twin.get_yscale())
                ylim = transform.inverted().transform(np.array(ylim_orig))
                if np.sign(np.diff(ylim_orig)) != np.sign(np.diff(ylim)): # dunno why needed
                    ylim = ylim[::-1]
                twin.set_ylim(offset + scale*ylim) # extra bit comes after the forward transformation

    def _rowlabels(self, labels, **kwargs):
        """Assign row labels."""
        axs = [ax for ax in self._main_axes if ax._xrange[0]==0]
        if isinstance(labels, str): # common during testing
            labels = [labels]*len(axs)
        if len(labels)!=len(axs):
            raise ValueError(f'Got {len(labels)} labels, but there are {len(axs)} rows.')
        axs = [ax for _,ax in sorted(zip([ax._yrange[0] for ax in axs], axs))]
        for ax,label in zip(axs,labels):
            if label and not ax.rowlabel.get_text():
                # Create a CompositeTransform that converts coordinates to
                # universal dots, then back to axes
                label_to_ax = ax.yaxis.label.get_transform() + ax.transAxes.inverted()
                x, _ = label_to_ax.transform(ax.yaxis.label.get_position())
                ax.rowlabel.update({'x':x, 'text':f'{label.strip()}   ', 'visible':True, **kwargs})

    def _collabels(self, labels, **kwargs):
        """Assign column labels."""
        axs = [ax for ax in self._main_axes if ax._yrange[0]==0]
        if isinstance(labels, str):
            labels = [labels]*len(axs)
        if len(labels)!=len(axs):
            raise ValueError(f'Got {len(labels)} labels, but there are {len(axs)} columns.')
        axs = [ax for _,ax in sorted(zip([ax._xrange[0] for ax in axs],axs))]
        for ax,label in zip(axs,labels):
            if all(pax.visible() for pax in ax.toppanel):
                ax = ax.toppanel[0]
            # if label and not ax.collabel.get_text():
            #     ax.collabel.update({'text':label, **kwargs})
            if label and not ax.collabel.get_text():
                ax.collabel.update({'text':label, **kwargs})

    def _suptitle_setup(self, renderer=None, **kwargs):
        """Intelligently determine super title position."""
        # Seems that `~matplotlib.figure.Figure.draw` is called *more than once*,
        # and the title position are appropriately offset only during the
        # *later* calls! So need to run this every time.
        # Row label; need to update as axis tick labels are generated and axis
        # label is offset!
        for ax in self._main_axes:
            if all(pax.visible() for pax in ax.toppanel):
                ax = ax.toppanel[0]
            label_to_ax = ax.yaxis.label.get_transform() + ax.transAxes.inverted()
            x, _ = label_to_ax.transform(ax.yaxis.label.get_position())
            ax.rowlabel.update({'x':x, 'ha':'center', 'transform':ax.transAxes})

        # Super title
        # Determine x by the underlying gridspec structure, where main axes lie.
        left = self._subplots_kw['left']
        right = self._subplots_kw['right']
        width, height = self.get_size_inches()
        if self.leftpanel:
            left += (self._subplots_kw['lwidth'] + self._subplots_kw['lspace'])
        if self.rightpanel:
            right += (self._subplots_kw['rwidth'] + self._subplots_kw['rspace'])
        xpos = left/width + 0.5*(width - left - right)/width

        # Figure out which title on the top-row axes will be offset the most
        # See: https://matplotlib.org/_modules/matplotlib/axis.html#Axis.set_tick_params
        # TODO: Modify self.axes? Maybe add a axes_main and axes_panel
        # attribute or something? Because if have insets and other stuff
        # won't that fuck shit up?
        # WARNING: Have to use private API to figure out whether axis has
        # tick labels or not! And buried in docs (https://matplotlib.org/api/_as_gen/matplotlib.axis.Axis.set_ticklabels.html)
        # label1 refers to left or bottom, label2 to right or top!
        # WARNING: Found with xtickloc='both', xticklabelloc='top' or 'both',
        # got x axis ticks position 'unknown'!
        # NOTE: Default linespacing is 1.2; it has no get, only a setter.
        # See: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
        line = 0
        xline = 0
        xlabel = 0 # room for x-label
        title = None
        raxs = [ax for ax in self._main_axes if ax._yrange[0]==0]
        raxs = [pax for ax in raxs for pax in ax.toppanel if pax.visible()] or raxs # filter to just these if top panel is enabled
        for axm in raxs:
            axs = [ax for ax in (axm, axm._twinx_child, axm._twiny_child) if ax]
            if all(pax.visible() for pax in axm.leftpanel): # note that if there are any top panels in the top row, they will of course not have their own panel attributes
                axs.extend(axm.leftpanel)
            if all(pax.visible() for pax in axm.rightpanel):
                axs.extend(axm.rightpanel)
            for ax in axs:
                if (ax.title.get_text().strip() and not ax._title_inside) or \
                   (ax.abc.get_text().strip() and not ax._abc_inside) or \
                   ax.collabel.get_text().strip():
                    line = 1.2
                    title = ax.title
                pos = ax.xaxis.get_ticks_position()
                if pos in ('top', 'unknown', 'default'):
                    label_top = pos!='default' and ax.xaxis.label.get_text().strip()
                    ticklabels_top = ax.xaxis._major_tick_kw['label2On']
                    if ticklabels_top:
                        if label_top: # offset above *tick labels* and *axis label*
                            # xline = max((xline, 2.8))
                            xline = max((xline, 3.2)) # empirically best
                        else:
                            xline = max((xline, 1.6))
                    xlabel = xline*(ax.xaxis.label.get_size()/72)/height
                    if (label_top or ticklabels_top) and ax.title.get_text().strip():
                        title = ax.title
                        line = 1.2
        # Default to whatever last axes was
        if not title:
            title = raxs[0].title # always will be non-None

        # Get the transformed position
        if self._suptitle_transform:
            transform = self._suptitle_transform
        else:
            transform = title.axes._title_pos_transform + self.transFigure.inverted()
        ypos = transform.transform(title.axes._title_pos_init)[1]
        line = line*(title.get_size()/72)/height + xlabel
        ypos = ypos + line
        # Update settings
        self._suptitle.update({
            'position': (xpos, ypos), 'transform': self.transFigure,
            'ha':'center', 'va':'bottom',
            **kwargs})

    def _auto_smart_tight_layout(self, renderer=None):
        """Conditionally call `~Figure.smart_tight_layout`."""
        # Only proceed if this wasn't already done, or user did not want
        # the figure to have tight boundaries.
        if not self._smart_tight_init or not (\
            self._smart_tight_outer or
            self._smart_tight_subplot or
            self._smart_tight_panel):
            return
        # Bail if we find cartopy axes
        # 1) If you used set_bounds to zoom into part of a cartopy projection,
        # this can erroneously identify invisible edges of map as being part of boundary
        # 2) If you have gridliner text labels, matplotlib won't detect them.
        if any(ax._cartopy_gridliner_on or ax._cartopy_limited_extent for ax in self._main_axes):
            warnings.warn('Tight subplots do not work with cartopy gridline labels or when zooming into a projection. Use tight=False in your call to subplots().')
            self._smart_tight_init = False
            return
        # Proceed
        if not self._silent:
            print('Adjusting layout.')
        self.smart_tight_layout(renderer)

    def _panel_tight_layout(self, side, panels, spans, renderer):
        """From list of panels and the axes coordinates spanned by the
        main axes, figure out the necessary new 'spacing' to apply to the inner
        gridspec object."""
        # Initial stuff
        pad = self._smart_panelpad
        if all(not panel for panel in panels): # exists, but may be invisible
            return
        elif not all(panel for panel in panels):
            raise ValueError('Either all or no axes in this row should have panels.')
        # Iterate through panels, get maximum necessary spacing; assign
        # equal spacing to all so they are still aligned
        space = []
        gspecs = [panel._panels_main_gridspec for panel in panels]
        if side in 'lr':
            ratios = [gs.get_width_ratios() for gs in gspecs]
        else:
            ratios = [gs.get_height_ratios() for gs in gspecs]
        for panel,span in zip(panels,spans):
            if not panel.visible() and not panel._colorbar_child:
                continue
            bbox = (panel._colorbar_child or panel).get_tightbbox(renderer)
            if side in 'lr':
                pspan = _intervalx_errfix(panel, bbox)
                # pspan = bbox.intervalx
            else:
                pspan = _intervaly_errfix(panel, bbox)
                # pspan = bbox.intervaly
            if side in 'tr':
                space.append((pspan[0] - span[1])/self.dpi)
            else:
                space.append((span[0] - pspan[1])/self.dpi)
        # Finally update the ratios, and fix
        # NOTE: Should just be able to modify mutable width and height ratios
        # list, do not need setter! The getter does not return a copy or anything.
        if not space:
            warnings.warn('All panels in this row or column are invisible. That is really weird.')
            return
        orig = [sum(iratios) for iratios in ratios]
        for panel,iratios in zip(panels,ratios):
            idx = (-2 if side in 'br' else 1)
            if panel._flush: # assume user means, always want panels touching, in spite of ticks, etc.
                ispace = 0
            else:
                ispace = max([0, iratios[idx] - min(space) + pad])
            iratios[idx] = ispace
        # Add back the difference due to changing spaces to the wratios and hratios
        # NOTE: Do not do this because will mess up implied aspect ratio maybe?
        # Not sure why this is not necessary.
        # for iorig,iratios in zip(orig,ratios):
        #     idx = (-3 if side in 'br' else 2)
        #     iratios[idx] -= (iorig - sum(iratios))

    def smart_tight_layout(self, renderer=None):
        """
        This is called automatically when `~Figure.draw` or
        `~Figure.savefig` are called,
        unless the user set `tight` to ``False`` in their original
        call to `subplots`.

        Conforms figure edges to tight bounding box around content, without
        screwing up subplot aspect ratios, empty spaces, panel sizes, and such.
        Also fixes inner spaces, so that e.g. axis tick labels and axis labels
        do not overlap with other axes (which was even harder to do and is
        insanely awesome).

        Parameters
        ----------
        renderer : None or `~matplotlib.backend_bases.RendererBase`, optional
            The backend renderer.
        """
        #----------------------------------------------------------------------#
        # Tick fudge factor
        #----------------------------------------------------------------------#
        # NOTE: Matplotlib majorTicks, minorTicks attributes contain *un-updated*
        # lists of ticks. To internally update the list and get actual ticks
        # in light of axis locator, use get_major_ticks and get_minor_ticks.
        # WARNING: Matplotlib seems to always draw "tight" bounding box around
        # tick marks, whether or not there are actually tick marks there! Must have
        # to do with get_tick_padding command; matplotlib tight layout just
        # queries that, which checks .majorTicks and .minorTicks, doesn't use getters.
        # WARNING: Could not figure out consistent recipe for calculating tick
        # pad errors in the seemingly infinite number of possible combos of
        # major and minor tick lengths, sides, label sides, and axis label sides.
        # New recommendation is to just use the 'wflush' and 'hflush' overrides
        # if user needs axes directly flush against each other, and in other
        # scenarios error won't matter much -- user can fudge subplotpad if necessary.
        # for ax in self._main_axes:
        #     iaxs = [ax, ax.leftpanel, ax.rightpanel, ax.bottompanel, ax.toppanel]
        #     for iax in iaxs:
        #         if not isinstance(iax, axes.XYAxes) or not iax.visible():
        #             continue
        #         # Store error in *points*, because we use it to adjust bounding
        #         # box span, whose default units are in dots.
        #         # Left right ticks
        #         if not iax.yaxis.label.get_text().strip():
        #             side = None
        #         else:
        #             side = iax.yaxis.label_position
        #         # *Always* pad minor tick error, for some reason
        #         pad = (0,0)
        #         ticks = [*iax.yaxis.minorTicks] # copy
        #         ticks_after = iax.yaxis.get_minor_ticks()
        #         if ticks:
        #             lbool, rbool = 1, 1
        #             if ticks_after: # only fix the side where ticks not actually present
        #                 lbool = int(not ticks_after[0].tick1On)
        #                 rbool = int(not ticks_after[0].tick2On)
        #             pad = 1 + ticks[0].get_tick_padding()
        #             pad = (lbool*pad*(side != 'left'), rbool*pad*(side != 'right'))
        #         iax._ytick_pad_error += np.array(pad)*self.dpi/72 # is left, right tuple
        #         # Pad major tick error only if ticks present, and label side
        #         # pad = (0,0)
        #         # ticks = iax.yaxis.majorTicks
        #         # if ticks:
        #         #     pad = ticks[0].get_tick_padding()
        #         #     pad = (pad*ticks[0].tick1On*(side != 'left'),
        #         #            pad*ticks[0].tick2On*(side != 'right'))
        #         # iax._ytick_pad_error += np.array(pad)*self.dpi/72 # is left, right tuple

        #----------------------------------------------------------------------#
        # Put tight box *around* figure
        #----------------------------------------------------------------------#
        subplots_kw = self._subplots_kw
        if self._smart_tight_outer:
            # Get bounding box that encompasses *all artists*, compare to bounding
            # box used for saving *figure*
            pad = self._smart_outerpad
            obbox = self.bbox_inches # original bbox
            if not renderer: # cannot use the below on figure save! figure becomes a special FigurePDF class or something
                renderer = self.canvas.get_renderer()
            bbox = self.get_tightbbox(renderer)
            ox, oy = obbox.intervalx, obbox.intervaly
            x, y = bbox.intervalx, bbox.intervaly
            # Apply new kwargs
            if np.any(np.isnan(x) | np.isnan(y)):
                warnings.warn('Bounding box has NaNs, cannot get outer tight layout.')
            else:
                left, bottom, right, top = x[0], y[0], ox[1]-x[1], oy[1]-y[1] # want *deltas*
                left   = subplots_kw['left'] - left + pad
                right  = subplots_kw['right'] - right + pad
                bottom = subplots_kw['bottom'] - bottom + pad
                top    = subplots_kw['top'] - top + pad + self._extra_pad
                subplots_kw.update({'left':left, 'right':right, 'bottom':bottom, 'top':top})

        #----------------------------------------------------------------------#
        # Prevent overlapping axis tick labels and whatnot *within* figure
        #----------------------------------------------------------------------#
        if self._smart_tight_subplot:
            # First for the main axes, then add in the panels
            xspans, yspans, xrange, yrange = _ax_props(self._main_axes, renderer)
            for side,panel in zip('lbr', (self.leftpanel, self.bottompanel, self.rightpanel)):
                # Only test panels on *inside* of subplot region
                # If have overlap panels due to content from outer stacked panel
                # your shit's fucked up. Anyway testing for that is more complicated
                # in the for rspace, for cspace loops below.
                # TODO: Check this for different orders.
                if not panel:
                    continue
                if panel._order=='C':
                    if side=='b':
                        panel = panel[:panel._n] # _n is number of columns
                    elif side=='l':
                        panel = panel[panel._n-1::panel._n]
                    else:
                        panel = panel[0::panel._n]
                else:
                    if side=='b':
                        panel = panel[0::panel._n] # _n is number of rows
                    elif side=='l':
                        panel = panel[-panel._n:]
                    else:
                        panel = panel[:panel._n]
                xs, ys, xr, yr = _ax_props(panel, renderer)
                xspans, yspans = np.vstack((xspans, xs)), np.vstack((yspans, ys))
                xrange, yrange = np.vstack((xrange, xr)), np.vstack((yrange, yr))
            xrange[:,1] += 1 # now, ids are not endpoint inclusive
            yrange[:,1] += 1
            # Construct "wspace" and "hspace" including panel spaces
            wspace_orig, hspace_orig = subplots_kw['wspace'], subplots_kw['hspace'] # originals
            if self.leftpanel:
                wspace_orig = [subplots_kw['lspace'], *wspace_orig]
            if self.rightpanel:
                wspace_orig = [*wspace_orig, subplots_kw['rspace']]
            if self.bottompanel:
                hspace_orig = [*hspace_orig, subplots_kw['bspace']]
            # Find groups of axes with touching left/right sides
            # We generate a list (xgroups), each element corresponding to a
            # column of *space*, containing lists of dictionaries with 'l' and
            # 'r' keys, describing groups of "touching" axes in that column
            xgroups = []
            ncols = self._main_axes[0]._ncols # includes panels?
            nrows = self._main_axes[0]._nrows
            for cspace in range(1, ncols): # spaces between columns
                groups = []
                for row in range(nrows):
                    filt = (yrange[:,0] <= row) & (row < yrange[:,1])
                    if sum(filt)<=1: # no interface here
                        continue
                    right, = np.where(filt & (xrange[:,0] == cspace))
                    left,  = np.where(filt & (xrange[:,1] == cspace))
                    if not (left.size==1 and right.size==1):
                        continue # normally both zero, but one can be non-zero e.g. if have left and bottom panel, empty space in corner
                    left, right = left[0], right[0]
                    added = False
                    for group in groups:
                        if left in group['l'] or right in group['r']:
                            group['l'].update((left,))
                            group['r'].update((right,))
                            added = True
                            break
                    if not added:
                        groups.append({'l':{left}, 'r':{right}}) # form new group
                xgroups.append(groups)
            # Find groups of axes with touching bottom/top sides
            ygroups = []
            for rspace in range(1, nrows): # spaces between rows
                groups = []
                for col in range(ncols):
                    filt = (xrange[:,0] <= col) & (col < xrange[:,1])
                    top,    = np.where(filt & (yrange[:,0] == rspace)) # gridspec indices go top-to-bottom, left-to-right
                    bottom, = np.where(filt & (yrange[:,1] == rspace))
                    if not (bottom.size==1 and top.size==1):
                        continue
                    bottom, top = bottom[0], top[0]
                    added = False
                    for group in groups:
                        if bottom in group['b'] or top in group['t']:
                            group['b'].update((bottom,))
                            group['t'].update((top,))
                            added = True
                            break
                    if not added:
                        groups.append({'b':{bottom}, 't':{top}}) # form new group
                ygroups.append(groups)
            # Correct wspace and hspace
            pad = self._smart_subplotpad
            if self._subplot_wflush:
                wspace = [0]*len(wspace_orig)
            else:
                wspace = []
                for space_orig, groups in zip(wspace_orig, xgroups):
                    if not groups:
                        wspace.append(space_orig)
                        continue
                    seps = [] # will use *maximum* necessary separation
                    for group in groups: # groups with touching edges
                        left = max(xspans[idx,1] for idx in group['l']) # axes on left side of column
                        right = min(xspans[idx,0] for idx in group['r']) # axes on right side of column
                        seps.append((right - left)/self.dpi)
                    wspace.append(max((0, space_orig - min(seps) + pad)))
            if self._subplot_hflush:
                hspace = [0]*len(hspace_orig)
            else:
                hspace = []
                for space_orig, groups in zip(hspace_orig, ygroups):
                    if not groups:
                        hspace.append(space_orig)
                        continue
                    seps = [] # will use *maximum* necessary separation
                    for group in groups: # groups with touching edges
                        bottom = min(yspans[idx,0] for idx in group['b'])
                        top = max(yspans[idx,1] for idx in group['t'])
                        seps.append((bottom - top)/self.dpi)
                    hspace.append(max((0, space_orig - min(seps) + pad)))
            # If had panels, need to pull out args
            lspace, rspace, bspace = 0, 0, 0 # does not matter if no panels
            if self.leftpanel:
                lspace, *wspace = wspace
            if self.rightpanel:
                *wspace, rspace = wspace
            if self.bottompanel:
                *hspace, bspace = hspace
            subplots_kw.update({'wspace':wspace, 'hspace':hspace,
                'lspace':lspace, 'rspace':rspace, 'bspace':bspace})
            # Update main gridspec to reflect new panel sizes
            self._main_gridspec.update()

        #----------------------------------------------------------------------#
        # The same, but for spaces between *axes panels*
        #----------------------------------------------------------------------#
        if self._smart_tight_panel and any(ax._panels_main_gridspec for ax in self._main_axes): # i.e. any axes has non-EmptyPanel panels
            # Get bounding boxes right now, otherwise they are queried a bunch
            # of times which may slow things down unnecessarily (benchmark?)
            # WARNING: Seems the erroneous padding is only issue for
            # subplots in the *main* gridspec; subplots inside inner
            # gridspecs to fine. Maybe depends on version though?
            axs = self._main_axes
            pairs = [(ax,ax.get_tightbbox(renderer)) for ax in axs]
            xspans = [_intervalx_errfix(ax,bbox) for ax,bbox in pairs]
            yspans = [_intervaly_errfix(ax,bbox) for ax,bbox in pairs]
            # xspans = [bbox.intervalx for ax,bbox in pairs]
            # yspans = [bbox.intervaly for ax,bbox in pairs]
            # Bottom, top
            for row in range(nrows):
                pairs = [(ax.bottompanel[0], yspan) for ax,yspan in zip(axs,yspans)
                        if ax._yrange[1]==row]
                if pairs:
                    self._panel_tight_layout('b', *zip(*pairs), renderer)
                pairs = [(ax.toppanel[-1], yspan) for ax,yspan in zip(axs,yspans)
                        if ax._yrange[0]==row]
                if pairs:
                    self._panel_tight_layout('t', *zip(*pairs), renderer)
            # Left, right
            for col in range(ncols):
                pairs = [(ax.leftpanel[-1], xspan) for ax,xspan in zip(axs,xspans)
                        if ax._xrange[0]==col]
                if pairs:
                    self._panel_tight_layout('l', *zip(*pairs), renderer)
                pairs = [(ax.rightpanel[0], xspan) for ax,xspan in zip(axs,xspans)
                        if ax._xrange[1]==col]
                if pairs:
                    self._panel_tight_layout('r', *zip(*pairs), renderer)
            # Reference axes
            axref = axs[self._ref_num - 1]
            gs = axref._panels_main_gridspec
            if not gs:
                wextra = 0
                hextra = 0
            else:
                wextra = sum(w for i,w in enumerate(gs.get_width_ratios())
                    if i!=2*int(bool(axref.leftpanel))) # main subplot index is position 2 of left panel is present, 0 if not
                hextra = sum(w for i,w in enumerate(gs.get_height_ratios())
                    if i!=2*int(bool(axref.toppanel))) # as for wextra
            subplots_kw.update({'wextra':wextra, 'hextra':hextra})
            # Update main gridspec to reflect new panel widths
            # TODO: Necessary?
            self._main_gridspec.update()

        #----------------------------------------------------------------------#
        # Finish
        #----------------------------------------------------------------------#
        # Parse arguments and update stuff
        # WARNING: Need to update gridspec twice. First to get the width
        # and height ratios including spaces, then after because gridspec
        # changes propagate only *down*, not up; will not be applied otherwise.
        figsize, gridspec_kw, _ = _subplots_kwargs(**subplots_kw)
        self._main_gridspec.update(**gridspec_kw)
        self._smart_tight_init = False
        self.set_size_inches(figsize)
        width, height = figsize
        # Update width and height ratios of axes panel subplotspec,
        # to reflect the new axes heights and widths (ratios are stored in physical units, inches)
        # NOTE: These will show us the *new* width and height ratios
        for ax in self._main_axes:
            # Outer ratios
            igs = ax._panels_main_gridspec
            if not igs:
                continue
            gs = self._main_gridspec
            idx = (2 if ax.toppanel else 0, 2 if ax.leftpanel else 0)
            wratios = gs.get_width_ratios()
            hratios = gs.get_height_ratios()
            # Ratios and 'extra' space for axes and its panels
            iwratios = igs.get_width_ratios() # gets my *custom* ratios
            ihratios = igs.get_height_ratios()
            ihextra = sum(h for i,h in enumerate(ihratios) if i!=idx[0])
            iwextra = sum(w for i,w in enumerate(iwratios) if i!=idx[1])
            # Calculate
            xrange = 2*np.array(ax._xrange) # endpoint inclusive
            yrange = 2*np.array(ax._yrange)
            fullwidth = sum(wratios[xrange[0]:xrange[1]+1]) # including panels!
            fullheight = sum(hratios[yrange[0]:yrange[1]+1])
            axwidth = fullwidth - iwextra
            axheight = fullheight - ihextra
            # Update axes panels
            iwratios[2 if ax.leftpanel else 0] = axwidth
            ihratios[2 if ax.toppanel else 0] = axheight
            igs.set_width_ratios(iwratios)
            igs.set_height_ratios(ihratios)
        # Final update, then update attributes, and we're good to go
        self._main_gridspec.update(**gridspec_kw)
        for ax in self._main_axes:
            width_new = np.diff(ax._position.intervalx)*width
            height_new = np.diff(ax._position.intervaly)*height
            ax.width, ax.height = width_new, height_new

    # @_counter
    def draw(self, renderer, *args, **kwargs):
        """Fix the "super title" position and automatically adjust the
        main gridspec, then call the parent `~matplotlib.figure.Figure.draw`
        method."""
        # Figure out if other titles are present, and if not bring suptitle
        # close to center. Also fix spanning labels (they need figure-relative
        # height coordinates, so must be updated when figure shape changes).
        self._lock_twins()
        self._suptitle_setup(renderer) # just applies the spacing
        self._auto_smart_tight_layout(renderer)
        for axis in self._span_labels:
            axis.axes._share_span_label(axis, final=True)
        out = super().draw(renderer, *args, **kwargs)

        # If rc settings have been changed, reset them after drawing is done
        # Usually means we have finished executing a notebook cell
        if not rc._init and self._rcreset:
            if not self._silent:
                print('Resetting rcparams.')
            rc.reset()
        return out

    def save(self, *args, **kwargs):
        """Alias for `~Figure.savefig`."""
        return self.savefig(*args, **kwargs)

    def savefig(self, filename, **kwargs):
        """
        Fix the "super title" position and automatically adjust the
        main gridspec, then call the parent `~matplotlib.figure.Figure.savefig`
        method.

        Parameters
        ----------
        filename : str
            The file name and path.
        **kwargs
            Passed to the matplotlib `~matplotlib.figure.Figure.savefig` method.
        """
        # Notes
        # * Gridspec object must be updated before figure is printed to
        #   screen in interactive environment; will fail to update after that.
        #   Seems to be glitch, should open thread on GitHub.
        # * To color axes patches, you may have to explicitly pass the
        #   transparent=False kwarg.
        #   Some kwarg translations, to pass to savefig
        if 'alpha' in kwargs:
            kwargs['transparent'] = not bool(kwargs.pop('alpha')) # 1 is non-transparent
        if 'color' in kwargs:
            kwargs['facecolor'] = kwargs.pop('color') # the color
            kwargs['transparent'] = False
        # Finally, save
        self._lock_twins()
        self._suptitle_setup() # just applies the spacing
        self._auto_smart_tight_layout()
        for axis in self._span_labels:
            axis.axes._share_span_label(axis, final=True)
        if not self._silent:
            print(f'Saving to "{filename}".')
        return super().savefig(os.path.expanduser(filename), **kwargs) # specify DPI for embedded raster objects

    def add_subplot_and_panels(self, subspec, which=None, *,
            hspace, wspace,
            bwidth, bvisible, bflush, bshare, bsep,
            twidth, tvisible, tflush, tshare, tsep,
            lwidth, lvisible, lflush, lshare, lsep,
            rwidth, rvisible, rflush, rshare, rsep,
            **kwargs):
        """
        Creates axes with optional "panels" along the sides of the axes. This
        is called by `subplots` -- you should not have to use this directly.

        Parameters
        ----------
        subspec : `~matplotlib.gridspec.SubplotSpec`
            The `~matplotlib.gridspec.SubplotSpec` instance onto which
            the main subplot and its panels are drawn.
        which : None or str, optional
            Whether to draw panels on the left, right, bottom, or top
            sides. Should be a string containing any of the characters
            ``'l'``, ``'r'``, ``'b'``, or ``'t'`` in any order.
            Default is ``'r'``.
        stack, lstack, rstack, bstack, tstack : int, optional
            The number of optional "stacked" panels on the left, right, bottom,
            and top sides, respectively. The default is ``1``. Use `stack` to
            set for all sides at once.
        sep, lsep, rsep, bsep, tsep : None, or float or str or list thereof, optional
            The separation between stacked panels. If float, units are inches.
            If string, units are interpreted by `~proplot.utils.units`. Ignored
            if the respecitve `stack` keyword arg is None. Use `sep` to set
            for all sides at once.
        share, lshare, rshare, bshare, tshare : bool, optional
            Whether to enable axis sharing between *y* axis of main subplot
            and left, right panels, and between *x* axis of main subplot and
            bottom, top panels. Use `share` to set for all
            sides at once.
        visible, lvisible, rvisible, bvisible, tvisible : bool, optional
            Used internally. Whether to make the left, right, bottom, or top
            panel invisible. This helps auto-align rows and columns of subplots
            when they don't all have panels on the same sides.
        width, lwidth, rwidth, bwidth, twidth : float or str or list thereof, optional
            Width of left, right, bottom, and top panels, respectively.
            If float or str, widths are same for all panels in the stack. If
            list thereof, specifies widths of each panel in the stack.
            If float, units are inches. If string, units are interpreted
            by `~proplot.utils.units`. Use `width` to set for
            all sides at once.
        wspace, hspace : float or str or list thereof, optional
            Empty space between the main subplot and the panel.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        flush, lflush, rflush, bflush, tflush : bool, optional
            Whether *inner* stacked panel should always be flush against the
            subplot, and *stacked* panels flush against each other.
            This overrides the `~Figure.smart_tight_layout` automatic spacing,
            and potentially overrides the `wspace`, `hspace`, and `lsep`,
            `rsep`, `bsep`, and `tsep` manual spacing. The default is ``False``.
            Use `flush` to set for all sides at once.
        """
        # Which panels
        # TODO: Allow visibility for some panels in stack, but not others?
        # This starts getting way too complicated.
        # NOTE: Axis sharing setup is handled after-the-fact in the `subplots`
        # function using `~proplot.axes.BaseAxes._sharex_setup` and
        # `~proplot.axes.BaseAxes._sharey_setup`.
        which = _default(which, 'r')
        if re.sub('[lrbt]', '', which): # i.e. other characters are present
            raise ValueError(f'Whichpanels argument can contain characters l (left), r (right), b (bottom), or t (top), instead got "{whichpanels}".')

        # Fix wspace/hspace in inches, using the Bbox from get_postition
        # on the subspec object to determine physical width of axes to be created
        # Fix widths and heights
        bbox = subspec.get_position(self) # valid since axes not drawn yet
        figwidth, figheight = self.get_size_inches()
        boxheight = np.diff(bbox.intervaly)[0]*figheight
        boxwidth = np.diff(bbox.intervalx)[0]*figwidth
        height = boxheight - sum(hspace)
        width = boxwidth - sum(wspace)
        # hspace = hspace/(height/nrows)
        # wspace = wspace/(width/ncols)

        # Figure out hratios/wratios
        # Will enforce (main_width + panel_width)/total_width = 1
        wratios = [width]
        if 'l' in which:
            lextra = sum(lwidth) + sum(lsep)
            wratios = [lextra, wratios[0] - lextra, *wratios[1:]]
        if 'r' in which:
            rextra = sum(rwidth) + sum(rsep)
            wratios = [*wratios[:-1], wratios[-1] - rextra, rextra]
        hratios = [height]
        if 't' in which:
            textra = sum(twidth) + sum(tsep)
            hratios = [textra, hratios[0] - textra, *hratios[1:]]
        if 'b' in which:
            bextra = sum(bwidth) + sum(bsep)
            hratios = [*hratios[:-1], hratios[-1] - bextra, bextra]
        if any(w < 0 for w in wratios):
            raise ValueError(f'Left and/or right panel widths too large for available space {width}in.')
        if any(h < 0 for h in hratios):
            raise ValueError(f'Top and/or bottom panel widths too large for available space {height}in.')

        # Make subplotspec
        panels = {}
        nrows = 1 + len(re.sub('[^bt]', '', which))
        ncols = 1 + len(re.sub('[^lr]', '', which))
        rows = np.array([w for w in ('t', 'c', 'b') if w in ('c', *which)])
        cols = np.array([w for w in ('l', 'c', 'r') if w in ('c', *which)])
        gs = gridspec.FlexibleGridSpecFromSubplotSpec(subplot_spec=subspec,
            nrows=nrows, ncols=ncols,
            wspace=wspace, hspace=hspace,
            width_ratios=wratios, height_ratios=hratios)
        # Draw main axes
        r, = np.where('c'==rows)
        c, = np.where('c'==cols)
        ax = self.add_subplot(gs[r[0],c[0]], **kwargs)
        ax._panels_main_gridspec = gs
        # Draw panels
        # TODO: Is this kwargs stuff necessary? What keyword args can user even supply?
        kwargs = {key:value for key,value in kwargs.items() if key not in ('number', 'projection')}
        for side,width,sep,flush,share,visible, in zip('lrbt',
                (lwidth, rwidth, bwidth, twidth),
                (lsep, rsep, bsep, tsep),
                (lflush, rflush, bflush, tflush),
                (lshare, rshare, bshare, tshare),
                (lvisible, rvisible, bvisible, tvisible),
                ):
            if side not in which:
                continue
            # Settings
            stack = len(width)
            if side in 'lr':
                r, = np.where('c'==rows)
                c, = np.where(side==cols)
                wspace, hspace = sep, []
                wratios, hratios = width, [1]
                nrows, ncols = 1, stack
            else:
                r, = np.where(side==rows)
                c, = np.where('c'==cols)
                wspace, hspace = [], sep
                wratios, hratios = [1], width
                nrows, ncols = stack, 1
            # Make gridspec and draw panels
            paxs = []
            name = {'b':'bottom', 't':'top', 'l':'left', 'r':'right'}[side]
            igs = gridspec.FlexibleGridSpecFromSubplotSpec(
                subplot_spec=gs[r[0],c[0]],
                nrows=nrows, ncols=ncols,
                wspace=wspace, hspace=hspace,
                width_ratios=wratios, height_ratios=hratios)
            for i in range(stack):
                pax = self.add_subplot(igs[i], side=name,
                    share=(visible and share), flush=flush, visible=visible,
                    parent=ax, projection='panel', **kwargs)
                pax._panels_main_gridspec = gs
                pax._panels_stack_gridspec = igs
                paxs += [pax]
            setattr(ax, name + 'panel', axes_list(paxs))
        # Set up axis sharing
        # * Sharex and sharey methods should be called on the 'child' axes,
        #   with labelling preferred on the axes in the first argument.
        # * This block must come *after* the above; sharex_setup and
        #   sharey_setup will detect invisible panels and not share with them.
        # First the x-axis
        bottom = None
        if 'b' in which and bvisible and bshare:
            bottom = ax.bottompanel[-1]
            for iax in (ax, *ax.bottompanel[:-1]):
                iax._sharex_setup(bottom, 3) # parent is *bottom-most* panel
        if 't' in which and tvisible and tshare:
            bottom = bottom or ax
            for iax in ax.toppanel:
                iax._sharex_setup(bottom, 3)
        # Now the y-axis
        left = None
        if 'l' in which and lvisible and lshare:
            left = ax.leftpanel[0]
            for iax in (*ax.leftpanel[1:], ax):
                iax._sharey_setup(left, 3) # parent is *bottom-most* panel
        if 'r' in which and rvisible and rshare:
            left = left or ax
            for iax in ax.rightpanel:
                iax._sharey_setup(left, 3)
        return ax

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
# Helper functions
def _panels_kwargs(panels, colorbars, legends, panels_kw, colorbars_kw=None, legends_kw=None,
        figure=False, ncols=None, nrows=None):
    """Returns standardized keyword args for axes and figure panels."""
    # Get which panels
    kwout = {}
    translate = {'bottom':'b', 'top':'t', 'left':'l', 'right':'r'}
    panels = translate.get(panels, panels)
    legends = translate.get(legends, legends)
    colorbars = translate.get(colorbars, colorbars)
    allpanels = panels + legends + colorbars
    allsides = 'lrb' if figure else 'lrbt'
    if len({*allpanels}) != len(allpanels):
        raise ValueError('You requested the same side for a panel, colorbar, and/or legend.')
    # Fill non-panels with empty args, copy over extra args to potentially
    # raise errors down the line.
    # NOTE: This is done to require keyword-only arguments to _subplots_kwargs
    # and add_subplot_and_panels, so we don't have to duplicate the rc[] stuff.
    if figure:
        names = ('width', 'sep', 'flush', 'share', 'space', 'span')
    else:
        names = ('width', 'sep', 'flush', 'share', 'visible')
    for side in {*allsides} - {*allpanels}: # add in specific ones
        for name in names:
            kwout[side + name] = None
    regex = re.compile(f'^[tlrb]?({"|".join(names)})$')
    if colorbars_kw is None:
        colorbars_kw = panels_kw
    if legends_kw is None:
        legends_kw = panels_kw
    for kwargs in (panels_kw, colorbars_kw, legends_kw):
        for key,value in kwargs.items():
            if not regex.match(key):
                kwout[key] = value
    # Helper function
    def get(side, name, defaults):
        """Get panel properties."""
        value = None
        if not isinstance(defaults, tuple):
            defaults = 3*(defaults,)
        for check,kwargs,default in zip((panels, colorbars, legends), (panels_kw, colorbars_kw, legends_kw), defaults):
            if side not in check:
                continue
            if isinstance(default, str):
                default = rc['subplot.' + default]
            return _default(kwargs.get(name, None), kwargs.get(side + name, None), default)

    # Get panel widths, account for stacked panels
    # NOTE: Accounts for stacked panels
    for side in allpanels:
        # Simple props
        stack = get(side, 'stack', 1)
        share = get(side, 'share', (True,False,True))
        kwout[side + 'share'] = share
        # Widths
        width = get(side, 'width', ('panelwidth', 'cbarwidth', 'legwidth'))
        width = np.atleast_1d(units(width))
        if len(width)==1:
            width = np.repeat(width, (stack,))
        if len(width)!=stack:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(width)} widths.')
        kwout[side + 'width'] = width
        # Panel separation
        # TODO: Figure out how to specify, e.g., a side panel and a side
        # colorbar together.
        flush = get(side, 'flush', False)
        if stack==1 or flush:
            sep = 0
        else:
            default = 'nolabspace' if share else 'ylabspace' if side in 'lr' else 'xlabspace'
            sep = get(side, 'sep', (default,default,default))
        sep = np.atleast_1d(units(sep))
        if len(sep)==1:
            sep = np.repeat(sep, (stack-1,))
        if len(sep)!=stack-1:
            raise ValueError(f'For side "{side}", have {stack} stacked panels, but got {len(sep)} separations.')
        kwout[side + 'sep'] = sep
        kwout[side + 'flush'] = flush

    if figure:
        for side in allpanels:
            # Space between panels and main subplots
            space = get(side, 'space', 'xlabspace' if side=='b' else 'ylabspace' if side=='l' else 'nolabspace')
            kwout[side + 'space'] = units(space)
            # Spanning of panels along subplot rows and columns
            nmax = ncols if side=='b' else nrows
            span = get(side, 'span', False)
            if np.iterable(span):
                ipanels = [*span]
                if len(span)!=nmax:
                    raise ValueError(f'Expected {nmax} {side}panel entries, got {len(span)}.')
            elif span:
                ipanels = [1]*nmax
            else:
                ipanels = [*range(1,nmax+1)]
            kwout[side + 'span'] = ipanels
    else:
        # Panel visibility, toggling
        kwout['which'] = panels + colorbars
        for side in allpanels:
            kwout[side + 'visible'] = get(side, 'visible', True)
        # Space between panels and parent subplot
        for name,axis,sides in zip(('hspace','wspace'), ('x','y'), ('tb','lr')):
            n = len([side for side in sides if side in allpanels])
            if not n:
                space = []
            else:
                big = axis + 'labspace'
                space = n*[0]
                for i,side in enumerate(sides):
                    if side not in allpanels:
                        continue
                    if kwout[side + 'flush']:
                        default = 0
                    else:
                        default = ('panelspace', big if side in 'bl' else 'panelspace', 'panelspace')
                    space[-i] = get(side, 'space', default) # not this also tests for 'lhspace', etc.; nonsense, but it's harmless
            kwout[name] = units(space)

    # Return the dictionary
    return kwout

def _subplots_kwargs(nrows, ncols, aspect, ref, *, # ref is the reference axes used for scaling things
    # Subplots
    # Empty args are filled with rc settings by main body of subplots()
    wextra, hextra,
    width,  height, axwidth, axheight,
    hspace, wspace, hratios, wratios,
    left,   bottom, right,   top,
    # Panels
    # Empty args are filled with rc settings by _panels_kwargs (helper function
    # meant to make figure panel, axes panel keyword args consistent)
    bspan, bwidth, bspace, bsep, bflush, bshare,
    lspan, lwidth, lspace, lsep, lflush, lshare,
    rspan, rwidth, rspace, rsep, rflush, rshare,
    ):
    """Handle complex keyword args and aliases thereof, and determine figure
    sizes such that aspect ratio of first subplot is preserved. Note
    bwidth, bspace, etc. must be supplied or will get error, but this should
    be taken care of by _parse_panels."""
    # Necessary arguments to reconstruct this grid, with defaults filled in
    subplots_kw = {'nrows': nrows, 'ncols': ncols, 'aspect': aspect, 'ref': ref,
        'width': width, 'height': height, 'axwidth': axwidth, 'axheight': axheight,
        'wextra': wextra, 'hextra': hextra, 'hspace': hspace, 'wspace': wspace, 'hratios': hratios, 'wratios': wratios,
        'left': left, 'bottom': bottom, 'right': right, 'top': top,
        'bspan': bspan, 'lspan': lspan, 'rspan': rspan,
        'bwidth': bwidth, 'bsep': bsep, 'bspace': bspace, 'bflush': bflush, 'bshare': bshare, # separation between panels
        'rwidth': rwidth, 'rsep': rsep, 'rspace': rspace, 'rflush': rflush, 'rshare': rshare,
        'lwidth': lwidth, 'lsep': lsep, 'lspace': lspace, 'lflush': lflush, 'lshare': lshare,
        }

    # Determine figure size
    # If width and height are not fixed, will scale them to preserve aspect
    # ratio of the first plot
    auto_both = (width is None and height is None)
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    auto_neither = (width is not None and height is not None)
    bpanel_space = sum(bwidth) + sum(bsep) + bspace if bspan else 0
    rpanel_space = sum(rwidth) + sum(rsep) + rspace if rspan else 0
    lpanel_space = sum(lwidth) + sum(lsep) + lspace if lspan else 0
    if np.iterable(aspect):
        aspect = aspect[0]/aspect[1]
    # Determine average axes widths/heights
    # Default behavior: axes average 2.0 inches wide
    # aspect = aspect*(hratios[0]/np.mean(hratios))
    if auto_width or auto_neither:
        axheights = height - top - bottom - sum(hspace) - bpanel_space
    if auto_height or auto_neither:
        axwidths = width - left - right - sum(wspace) - rpanel_space - lpanel_space
    # If both are auto (i.e. no figure dims specified), this is where
    # we set figure dims according to axwidth and axheight input
    if auto_both: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = units(rc['subplot.axwidth'])
        if axheight is not None:
            height = axheight*nrows + top + bottom + sum(hspace) + bpanel_space
            axheights = axheight*nrows
            auto_width = True
        if axwidth is not None:
            width = axwidth*ncols + left + right + sum(wspace) + rpanel_space + lpanel_space
            axwidths = axwidth*ncols
            auto_height = True
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
        figsize = (width, height) # again

    # Automatically scale fig dimensions
    # TODO: The reference axes shit fails for complex grids!
    # TODO: Maybe it should refer to a *position on the grid*, since super hard
    # to generalize for arbitrary axes occupying arbitrary boxes on the grid?
    href = (ref - 1) % ncols
    wref = (ref - 1) // ncols
    if auto_width:
        axheight = axheights*hratios[href]/sum(hratios)
        axwidth  = (axheight - hextra)*aspect + wextra # bigger w ratio, bigger width
        axwidths = axwidth*sum(wratios)/wratios[wref]
        width = axwidths + left + right + sum(wspace) + rpanel_space + lpanel_space
    elif auto_height:
        axwidth = axwidths*wratios[wref]/sum(wratios)
        axheight = (axwidth - wextra)/aspect + hextra
        axheights = axheight*sum(hratios)/hratios[href]
        height = axheights + top + bottom + sum(hspace) + bpanel_space
    if axwidths<0:
        raise ValueError(f"Not enough room for axes (would have width {axwidths}). Try using tight=False, increasing figure width, or decreasing 'left', 'right', or 'wspace' spaces.")
    if axheights<0:
        raise ValueError(f"Not enough room for axes (would have height {axheights}). Try using tight=False, increasing figure height, or decreasing 'top', 'bottom', or 'hspace' spaces.")
    # Make sure the 'ratios' and 'spaces' are in physical units (we cast the
    # former to physical units), easier then to add stuff as below
    wspace = [*wspace]
    hspace = [*hspace]
    wratios = [*(axwidths*wratios/sum(wratios))]
    hratios = [*(axheights*hratios/sum(hratios))]
    # Now add the outer panel considerations, will enforce panels with
    # constant widths.
    nrows += int(bool(bspan))
    ncols += int(bool(rspan)) + int(bool(lspan))
    if bspan:
        hratios = hratios + [sum(bwidth) + sum(bsep)]
        hspace  = hspace  + [bspace]
    if lspan:
        wratios = [sum(lwidth) + sum(lsep)] + wratios
        wspace  = [lspace] + wspace
    if rspan:
        wratios = wratios + [sum(rwidth) + sum(rsep)]
        wspace  = wspace  + [rspace]
    # Keyword args for gridspec class
    bottom = bottom/height
    left   = left/width
    top    = 1 - top/height
    right  = 1 - right/width
    gridspec_kw = {
        'nrows': nrows, 'ncols': ncols,
        'left': left, 'bottom': bottom, 'right': right, 'top': top, # so far no panels allowed here
        'wspace': wspace, 'hspace': hspace, 'width_ratios': wratios, 'height_ratios' : hratios,
        }
    return (width, height), gridspec_kw, subplots_kw

def _axes_dict(naxs, value, kw=False, default=None):
    """Build a dictionary that looks like {1:value1, 2:value2, ...} or
    {1:{key1:value1, ...}, 2:{key2:value2, ...}, ...} for storing
    standardized axes-specific properties or keyword args."""
    # First build up dictionary
    # 1) 'string' or {1:'string1', (2,3):'string2'}
    if not kw:
        if np.iterable(value) and not isinstance(value, (str,dict)):
            value = {num+1: item for num,item in enumerate(value)}
        elif not isinstance(value, dict):
            value = {range(1, naxs+1): value}
    # 2) {'prop':value} or {1:{'prop':value1}, (2,3):{'prop':value2}}
    else:
        nested = [isinstance(value,dict) for value in value.values()]
        if not any(nested): # any([]) == False
            value = {range(1, naxs+1): value.copy()}
        elif not all(nested):
            raise ValueError('Pass either of dictionary of key value pairs or a dictionary of dictionaries of key value pairs.')
    # Then *unfurl* keys that contain multiple axes numbers, i.e. are meant
    # to indicate properties for multiple axes at once
    kwargs = {}
    for nums,item in value.items():
        nums = np.atleast_1d(nums)
        for num in nums.flat:
            if not kw:
                kwargs[num] = item
            else:
                kwargs[num] = item.copy()
    # Fill with default values
    for num in range(1, naxs+1):
        if num not in kwargs:
            if kw:
                kwargs[num] = {}
            else:
                kwargs[num] = default
    # Verify numbers
    if {*range(1, naxs+1)} != {*kwargs.keys()}:
        raise ValueError(f'Have {naxs} axes, but {value} has properties for axes {", ".join(str(i+1) for i in sorted(kwargs.keys()))}.')
    return kwargs

# Primary function
def subplots(array=None, ncols=1, nrows=1,
        # Figure settings
        ref=1, # reference axes for fixing aspect ratio
        order='C', # allow calling with subplots(array)
        aspect=1, figsize=None,
        width=None,  height=None, axwidth=None, axheight=None, journal=None,
        wwidth=None, hwidth=None,
        axwidths=None, axheights=None,
        hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
        width_ratios=None, height_ratios=None,
        flush=False, wflush=None, hflush=None,
        left=None, bottom=None, right=None, top=None, # spaces around edge of main plotting area, in inches
        # Padding settings
        emptycols=None, emptyrows=None, # obsolete?
        tight=None, outertight=None, subplottight=None, paneltight=None,
        outerpad=None, panelpad=None, subplotpad=None,
        # Shared axes and spanning labels
        span=None,  spanx=1,  spany=1,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        share=None, sharex=3, sharey=3, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        # Panels and projections
        panel=None, panels=None, legend=None, legends=None, colorbar=None, colorbars=None,
        axpanel=None, axlegend=None, axcolorbar=None,
        axpanels=None, axlegends=None, axcolorbars=None,
        axpanel_kw=None, axcolorbar_kw=None, axlegend_kw=None,
        axpanels_kw=None, axcolorbars_kw=None, axlegends_kw=None,
        basemap=False, proj=None, projection=None, proj_kw=None, projection_kw=None,
        rcreset=True, silent=True, # arguments for figure instantiation
        **kwargs):
    """
    Analagous to `matplotlib.pyplot.subplots`. Create a figure with a single
    axes or arbitrary grids of axes (any of which can be map projections),
    and optional "panels" along axes or figure edges.

    The parameters are sorted into the following rough sections: subplot grid
    specifications, figure and subplot sizes, axis sharing,
    figure panels, axes panels, and map projections.

    Parameters
    ----------
    ncols, nrows : int, optional
        Number of columns, rows. Ignored if `array` is not ``None``.
        Use these arguments for simpler subplot grids.
    order : {'C', 'F'}, optional
        Whether subplots are numbered in column-major (``'C'``) or row-major
        (``'F'``) order. Analagous to `numpy.array` ordering. This controls
        the order axes appear in the `axs` list, and the order of a-b-c
        labelling when using `~proplot.axes.BaseAxes.format` with ``abc=True``.
    array : None or array-like of int, optional
        2-dimensional array specifying complex grid of subplots. Think of
        this array as a "picture" of your figure. For example, the array
        ``[[1, 1], [2, 3]]`` creates one long subplot in the top row, two
        smaller subplots in the bottom row.

        Integers must range from 1 to the number of plots. For example,
        ``[[1, 4]]`` is invalid.

        ``0`` indicates an empty space. For example, ``[[1, 1, 1], [2, 0, 3]]``
        creates one long subplot in the top row with two subplots in the bottom
        row separated by a space.
    emptyrows, emptycols : None or list of int, optional
        Row, column numbers (starting from 1) that you want empty. Generally
        this is used with `ncols` and `nrows`; with `array`, you can
        just use zeros.

    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    journal : None or str, optional
        Conform figure width (and height, if specified) to academic journal
        standards. See `~proplot.gridspec.journal_size` for details.
    width, height : float or str, optional
        The figure width and height. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units` -- for example,
        `width="10cm"` creates a 10cm wide figure.
    axwidth, axheight : float or str, optional
        Sets the average width, height of your axes. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.

        These arguments are convenient where you don't care about the figure
        dimensions and just want your axes to have enough "room".
    aspect : float or length-2 list of floats, optional
        The (average) axes aspect ratio, in numeric form (width divided by
        height) or as (width, height) tuple. If you do not provide
        the `hratios` or `wratios` keyword args, all axes will have
        identical aspect ratios.
    hratios, wratios
        Shorthands for `height_ratios`, `width_ratios`.
    hspace, wspace, height_ratios, width_ratios, left, right, top, bottom : optional
        Passed to `~proplot.gridspec.FlexibleGridSpecBase`.

    sharex, sharey, share : {3, 2, 1, 0}, optional
        The "axis sharing level" for the *x* axis, *y* axis, or both
        axes. Options are as follows:

            0. No axis sharing.
            1. Only draw *axis label* on the leftmost column (*y*) or
               bottommost row (*x*) of subplots. Axis tick labels
               still appear on every subplot.
            2. As in 1, but forces the axis limits to be identical. Axis
               tick labels still appear on every subplot.
            3. As in 2, but only show the *axis tick labels* on the
               leftmost column (*y*) or bottommost row (*x*) of subplots.

        This feature can considerably redundancy in your figure.
    spanx, spany, span : bool or {1, 0}, optional
        Whether to use "spanning" axis labels for the *x* axis, *y* axis, or
        both axes. When ``True`` or ``1``, the axis label for the leftmost
        (*y*) or bottommost (*x*) subplot is **centered** on that column or
        row -- i.e. it "spans" that column or row. For example, this means
        for a 3-row figure, you only need 1 y-axis label instead of 3.

        This feature can considerably redundancy in your figure.

        "Spanning" labels also integrate with "shared" axes. For example,
        for a 3-row, 3-column figure, with ``sharey>1`` and ``spany=1``,
        your figure will have **1** *y*-axis label instead of 9.

    proj, proj_kw
        Aliases for `projection`, `projection_kw`.
    projection : str or dict-like, optional
        The map projection name.

        If string, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `axpanels`.

        For example, consider a figure with 4 columns and 1 row.
        With ``projection={1:'mercator', (2,3):'hammer'}``,
        the leftmost subplot is a Mercator projection, the middle 2 are Hammer
        projections, and the rightmost is a normal x/y axes.
    projection_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class on instantiation.

        As with `axpanels_kw`, can be dict of properties (applies globally),
        or **dict of dicts** of properties (applies to specific axes, as
        with `axpanels`).
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.

        If boolean, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `axpanels`.

    panel : str, optional
        Specify which sides of the figure should have a "panel".
        This is great for creating global colorbars and legends.
        String should contain any of the characters ``'l'`` (left panel),
        ``'r'`` (right panel), or ``'b'`` (bottom panel). For example,
        ``'br'`` will draw a right and bottom panel.

        Panels are stored on the ``leftpanel``, ``rightpanel``,
        and ``bottompanel`` attributes on the figure instance.
        You can also use the aliases ``lpanel``, ``rpanel``, and ``bpanel``.
    panels : str, optional
        Same as `panel`, but default behavior is to assign a panel
        to *every* row or column of subplots. Individual panels can then
        be accessed with e.g. ``fig.leftpanel[0]``, ``fig.leftpanel[1]``.
    colorbar, legend, colorbars, legends
        Identical to `panel` and `panels`, except the *default* panel width is
        more appropriate for a colorbar or legend. The panel can then be
        "filled" with a colorbar or legend with e.g.
        ``fig.bottompanel.colorbar()`` or ``fig.bottompanel.legend()``.
    lspan, rspan, bspan : bool or list of int, optional
        Define how figure panels span rows and columns of subplots.
        Argument is interpreted as follows:

        * If ``True``, this is default behavior for `panel` keyword arg. Draws
          single panel spanning **all** columns/rows of subplots.
        * If ``False``, this is default behavior for `panels` keyword arg. Draws
          separate panels **for each** column/row of subplots.
        * If list of int, you can specify panels that span **contiguous**
          columns/rows of subplots. Usage is similar to the ``array`` argument.

          For example, for a figure with 3 columns, ``bspan=[1, 2, 2]``
          draws a panel on the bottom of the first column and spanning the
          bottom of the right 2 columns, and ``bspan=[0, 2, 2]``
          only draws a panel underneath the right 2 columns (as with
          `array`, the ``0`` indicates an empty space).
    lspace, rspace, bspace : float, optional
        Space between the "inner" edges of the left, right, and bottom
        panels and the edge of the main subplot grid. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.
    lwidth, rwidth, bwidth
        See `~Figure.add_subplot_and_panels`. Usage is identical.
    lstack, rstack, bstack
        See `~Figure.add_subplot_and_panels`. Usage is identical.
    lsep, rsep, bsep
        See `~Figure.add_subplot_and_panels`. Usage is identical.
    lflush, rflush, bflush
        As in `~Figure.add_subplot_and_panels`, but only controls whether
        *stacked* panels are flush against each other.
    lshare, rshare, bshare
        As in `~Figure.add_subplot_and_panels`, but only controls axis
        sharing between *stacked* panels.

    axpanel, axpanels : str or dict-like, optional
        Specify which axes should have "panels".
        Panels are stored on the ``leftpanel``, ``rightpanel``,
        ``bottompanel`` and ``toppanel`` attributes on the axes instance.
        You can also use the aliases ``lpanel``, ``rpanel``, ``bpanel``,
        or ``tpanel``.

        The argument is interpreted as follows:

            * If str, panels are drawn on the same side for all subplots.
              String should contain any of the characters ``'l'`` (left panel),
              ``'r'`` (right panel), ``'t'`` (top panel), or ``'b'`` (bottom panel).
              For example, ``'rt'`` will draw a right and top panel.
            * If dict-like, panels can be drawn on different sides for
              different subplots. For example, for a 4-subplot figure,
              ``axpanels={1:'r', (2,3):'l'}`` indicates that we want to
              draw a panel on the right side of subplot number 1, on the left
              side of subplots 2 and 3, and **no panel** on subplot 4.

    axcolorbar, axlegend, axcolorbars, axlegends
        Identical to `axpanel` or `axpanels`, except the *default* panel width is
        more appropriate for a colorbar or legend. The panel can then be
        "filled" with a colorbar or legend with e.g.
        ``ax.rightpanel.colorbar()`` or ``ax.rightpanel.legend()``.
    axpanel_kw, axcolorbar_kw, axlegend_kw, axcolorbars_kw, axlegends_kw
        Aliases for `axpanels_kw`.
    axpanels_kw : dict-like, optional
        Keyword args passed to `~Figure.add_subplot_and_panels`.
        Can be dict of properties (applies globally), or **dict of dicts** of
        properties (applies to specific axes, as with `axpanels`).

        For example, consider a 2-subplot figure with ``axpanels='l'``.
        With ``{'lwidth':1}``, both left panels will be 1 inch wide.
        With ``{1:{'lwidth':1}, 2:{'lwidth':0.5}}``, the left subplot
        panel will be 1 inch wide and the right subplot panel will be
        0.5 inches wide.

    Returns
    -------
    f : `Figure`
        The figure instance.
    axs : `axes_list`
        A special list of axes instances. See `axes_list`.

    Other parameters
    ----------------
    tight, outertight, subplottight, paneltight, outerpad, subplotpad, panelpad, flush, wflush, hflush, rcreset, silent
        Passed to `Figure`.
    **kwargs
        Passed to `~proplot.gridspec.FlexibleGridSpec`.

    See also
    --------
    `~proplot.axes.BaseAxes`, `~proplot.axes.BaseAxes.format`
    """
    # TODO: Make axpanel_kw and axcolorbar_kw apply to those separate keyword
    # args, optionally.
    # TODO: Generalize axes sharing for right y-axes and top x-axes. Enable a secondary
    # axes sharing mode where we *disable ticklabels and labels*, but *do not
    # use the builtin sharex/sharey API*, suitable for complex map projections.
    # For spanning axes labels, right now only detect **x labels on bottom**
    # and **ylabels on top**. Generalize for all subplot edges.
    #--------------------------------------------------------------------------#
    # Create blank figure
    # Will adjust dimensions and stuff later
    #--------------------------------------------------------------------------#
    fig = plt.figure(FigureClass=Figure, tight=tight,
        outertight=outertight, subplottight=subplottight, paneltight=paneltight,
        outerpad=outerpad,     subplotpad=subplotpad,     panelpad=panelpad,
        flush=flush, wflush=wflush, hflush=hflush, rcreset=rcreset, silent=silent,
        )

    #--------------------------------------------------------------------------#
    # Array setup
    #--------------------------------------------------------------------------#
    rc._getitem_mode = 0 # ensure still zero; might be non-zero if had error in 'with context' block
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        if order not in ('C','F'): # better error message
            raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
        array = array.reshape((nrows, ncols), order=order) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        if order not in ('C','F'): # better error message
            raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
        array = array[None,:] if order=='C' else array[:,None] # interpret as single row or column
    # Empty rows/columns feature
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    if emptycols:
        emptycols = np.atleast_1d(emptycols)
        for col in emptycols.flat:
            array[:,col-1] = 0
    if emptyrows:
        emptyrows = np.atleast_1d(emptyrows)
        for row in emptyrows.flat:
            array[row-1,:] = 0
    # Enforce rule
    nums = np.unique(array[array!=0])
    naxs = len(nums)
    if {*nums.flat} != {*range(1, naxs+1)}:
        raise ValueError('Axes numbers must span integers 1 to naxs (i.e. cannot skip over numbers).')
    nrows = array.shape[0]
    ncols = array.shape[1]

    #--------------------------------------------------------------------------#
    # Panels
    #--------------------------------------------------------------------------#
    # Parse outer panel kwargs, consolidate settings
    for ipanel in (panel,legend,colorbar):
        for side in (ipanel or ''):
            if side + 'span' not in kwargs:
                kwargs[side + 'span'] = True # spans all rows or columns
    panels    = _default(panel, panels, '')
    legends   = _default(legend, legends, '')
    colorbars = _default(colorbar, colorbars, '')
    # Get panel props
    kwargs = _panels_kwargs(panels, colorbars, legends, kwargs,
            figure=True, ncols=ncols, nrows=nrows)
    # Create dictionary of panel toggles and settings
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    # TODO: Allow separate settings for separate colorbar, legend, etc. panels
    axpanels    = _axes_dict(naxs, _default(axpanel, axpanels, ''), kw=False, default='')
    axcolorbars = _axes_dict(naxs, _default(axcolorbar, axcolorbars, ''), kw=False, default='')
    axlegends   = _axes_dict(naxs, _default(axlegend, axlegends, ''), kw=False, default='')
    axpanels_kw    = _axes_dict(naxs, _default(axpanel_kw, axpanels_kw, {}), kw=True)
    axcolorbars_kw = _axes_dict(naxs, _default(axcolorbar_kw, axcolorbars_kw, {}), kw=True)
    axlegends_kw   = _axes_dict(naxs, _default(axlegend_kw, axlegends_kw, {}), kw=True)
    # Get which panels
    for num in range(1,naxs+1):
        axpanels_kw[num] = _panels_kwargs(
            axpanels[num], axcolorbars[num], axlegends[num],
            axpanels_kw[num], axcolorbars_kw[num], axlegends_kw[num],
            figure=False)

    #--------------------------------------------------------------------------#
    # Shared and spanning axes
    #--------------------------------------------------------------------------#
    # Figure out rows and columns "spanned" by each axes in list, for
    # axis sharing and axis label spanning settings
    sharex = int(_default(share, sharex))
    sharey = int(_default(share, sharey))
    spanx  = _default(span, spanx)
    spany  = _default(span, spany)
    if sharex not in range(4) or sharey not in range(4):
        raise ValueError('Axis sharing options sharex/sharey can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels).')
    # Get some axes properties. Locations are sorted by axes id.
    axids = [np.where(array==i) for i in np.sort(np.unique(array)) if i>0] # 0 stands for empty
    idxoff = (0, int(bool(kwargs.get('lspan', None))))
    yrange = idxoff[0] + np.array([[y.min(), y.max()+1] for y,_ in axids]) # yrange is shared columns
    xrange = idxoff[1] + np.array([[x.min(), x.max()+1] for _,x in axids])

    # Shared axes: Generate list of base axes-dependent axes pairs
    # That is, find where the minimum-maximum gridspec extent in 'x' for a
    # given axes matches the minimum-maximum gridspec extent for a base axes.
    xgroups_base, xgroups, grouped = [], [], {*()}
    if sharex:
        for i in range(naxs): # axes now have pseudo-numbers from 0 to naxs-1
            matches       = (xrange[i,:]==xrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0] # gives ID number of matching_axes, from 0 to naxs-1
            if i not in grouped and matching_axes.size>1:
                # Find all axes that have the same gridspec 'x' extents
                xgroups += [matching_axes]
                # Get bottom-most axis with shared x; should be single number
                xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
            grouped.update(matching_axes) # bookkeeping; record ids that have been grouped already
    ygroups_base, ygroups, grouped = [], [], {*()}
    # Similar for *y* axis
    if sharey:
        for i in range(naxs):
            matches       = (yrange[i,:]==yrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups += [matching_axes]
                ygroups_base += [matching_axes[np.argmin(xrange[matching_axes,0])]] # left-most axis with shared y, for matching_axes
            grouped.update(matching_axes) # bookkeeping; record ids that have been grouped already

    # Make sure same panels are present on ***all rows and columns***
    # Example: If left panel was requested for top subplot in 2-row figure, will
    # allocate space for empty left panel in bottom subplot too, to keep
    # things aligned!
    # WARNING: This is very limited right now. Will just use the lowest number
    # in the dictionary (since np.where return is sorted, so idxs is sorted,
    # so filt_kw and side_kw are sorted, since python 3.6+ dicts preserve order)
    # WARNING: Used to test for maximum spacing, panel widths along row, but
    # now that have stacked panels this becomes just too complex. Note smart
    # tight layout does its own thing to ensure appropriate spacing.
    def sync(side, idxs):
        all_kw = {idx+1: axpanels_kw[idx+1] for idx in idxs}
        on_kw = {num:kw for num,kw in all_kw.items() if side in kw['which']}
        if not on_kw:
            return
        dict_ = (*on_kw.values(),)[0]
        idx = 0 if side in 'tl' else -1
        space_ = 'wspace' if side in 'lr' else 'hspace'
        sep, space, width = dict_[side + 'sep'], dict_[space_], dict_[side + 'width']
        for num,kw in all_kw.items():
            kw[side + 'sep'] = sep
            kw[side + 'width'] = width
            if num in on_kw:
                kw[space_][idx] = space[idx] # enforce, e.g. in case one left panel has different width compared to another!
                continue
            kw['which'] += side
            kw[side + 'visible'] = False
            if side in 'tl':
                kw[space_] = [space[idx], *kw[space_]]
            else:
                kw[space_] = [*kw[space_], space[idx]]
    for row in range(nrows):
        idxs, = np.where(yrange[:,0]==row) # note ranges are *slices*, not endpoint inclusive
        sync('t', idxs)
        idxs, = np.where((yrange[:,1]-1)==row) # note ranges are *slices*, not endpoint inclusive
        sync('b', idxs)
    for col in range(ncols):
        idxs, = np.where(xrange[:,0]==col) # note ranges are *slices*, not endpoint inclusive
        sync('l', idxs)
        idxs, = np.where((xrange[:,1]-1)==col) # note ranges are *slices*, not endpoint inclusive
        sync('r', idxs)

    #--------------------------------------------------------------------------#
    # Get basemap.Basemap or cartopy.CRS instances for map, and
    # override aspect ratio.
    #--------------------------------------------------------------------------#
    # NOTE: Cannot have mutable dict as default arg, because it changes the
    # "default" if user calls function more than once! Swap dicts for None.
    basemap = _axes_dict(naxs, basemap, kw=False, default=False)
    proj    = _axes_dict(naxs, _default(proj, projection, 'xy'), kw=False, default='xy')
    proj_kw = _axes_dict(naxs, _default(proj_kw, projection_kw, {}), kw=True)
    axes_kw = {num:{} for num in range(1, naxs+1)}  # stores add_subplot arguments
    for num,name in proj.items():
        # The default, my XYAxes projection
        if name=='xy':
            axes_kw[num]['projection'] = 'xy'
        # Builtin matplotlib polar axes, just use my overridden version
        elif name=='polar':
            axes_kw[num]['projection'] = 'newpolar'
            if num==ref:
                aspect = 1
        # Custom Basemap and Cartopy axes
        elif name:
            package = 'basemap' if basemap[num] else 'cartopy'
            instance, iaspect, kwproj = projs.Proj(name, basemap=basemap[num], **proj_kw[num])
            if num==ref:
                aspect = iaspect
            axes_kw[num].update({'projection':package, 'map_projection':instance})
            axes_kw[num].update(kwproj)
        else:
            raise ValueError('All projection names should be declared. Wut.')

    #--------------------------------------------------------------------------#
    # Figure architecture
    #--------------------------------------------------------------------------#
    # Figure and/or average axes dimensions
    if journal:
        figsize = journals(journal) # if user passed width=<string>, will use that journal size
    elif not figsize:
        figsize = (width, height)
    width, height = figsize
    width  = units(width)
    height = units(height)
    axwidth  = units(axwidth)
    axheight = units(axheight)
    # Gridspec defaults
    # NOTE: Ratios are scaled to take physical units in _subplots_kwargs, so
    # user can manually provide hspace and wspace in physical units.
    hspace = np.atleast_1d(units(_default(hspace,
        rc['subplot.titlespace'] + rc['subplot.innerspace'] if sharex==3
        else rc['subplot.xlabspace'] if sharex in (1,2) # space for tick labels and title
        else rc['subplot.titlespace'] + rc['subplot.xlabspace'])))
    wspace = np.atleast_1d(units(_default(wspace,
        rc['subplot.innerspace'] if sharey==3
        else rc['subplot.ylabspace'] - rc['subplot.titlespace'] if sharey in (1,2) # space for tick labels only
        else rc['subplot.ylabspace']
        )))
    wratios = np.atleast_1d(_default(width_ratios, wratios, 1))
    hratios = np.atleast_1d(_default(height_ratios, hratios, 1))
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    if len(wratios)==1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios)==1:
        hratios = np.repeat(hratios, (nrows,))
    left   = units(_default(left,   rc['subplot.ylabspace']))
    bottom = units(_default(bottom, rc['subplot.xlabspace']))
    right  = units(_default(right,  rc['subplot.nolabspace']))
    top    = units(_default(top,    rc['subplot.titlespace']))
    # Extra space necessary to account for fixing axes aspect ratios when
    # inner panels are present
    axkw = axpanels_kw[ref]
    wextra = sum(axkw['wspace']) + \
            sum(sum(_default(axkw[side + 'width'], [0])) for side in 'lr') + \
            sum(sum(_default(axkw[side + 'sep'], [0])) for side in 'lr')
    hextra = sum(axkw['hspace']) + \
            sum(sum(_default(axkw[side + 'width'], [0])) for side in 'tb') + \
            sum(sum(_default(axkw[side + 'sep'], [0])) for side in 'tb')
    # Parse arguments, fix dimensions in light of desired aspect ratio
    figsize, gridspec_kw, subplots_kw = _subplots_kwargs(nrows, ncols, aspect, ref,
            left=left, right=right, bottom=bottom, top=top,
            width=width, height=height, axwidth=axwidth, axheight=axheight,
            wratios=wratios, hratios=hratios, wspace=wspace, hspace=hspace, wextra=wextra, hextra=hextra,
            **kwargs)
    # Apply settings and add attributes
    gs = gridspec.FlexibleGridSpec(**gridspec_kw)
    fig.set_size_inches(figsize)
    fig._main_gridspec = gs
    fig._subplots_kw = subplots_kw

    #--------------------------------------------------------------------------#
    # Draw on figure
    #--------------------------------------------------------------------------#
    # Create axes
    axs = naxs*[None] # list of axes
    for idx in range(naxs):
        num = idx + 1
        ax_kw = axes_kw[num]
        if axpanels_kw[num]['which']: # non-empty
            axs[idx] = fig.add_subplot_and_panels(gs[slice(*yrange[idx,:]), slice(*xrange[idx,:])],
                    number=num, spanx=spanx, spany=spany, **ax_kw, **axpanels_kw[num])
        else:
            axs[idx] = fig.add_subplot(gs[slice(*yrange[idx,:]), slice(*xrange[idx,:])],
                    number=num, spanx=spanx, spany=spany, **ax_kw) # main axes can be a cartopy projection
    # Create outer panels
    for side in 'blr':
        # Draw panel from gridspec
        panels = subplots_kw[side + 'span']
        if not panels:
            continue
        paxs = []
        name = {'b':'bottom', 'l':'left', 'r':'right'}[side]
        for num in np.unique(panels).flat:
            # Get subspec, and settings
            # TODO: Allow for different width ratios and stuff
            if num==0:
                continue
            off = idxoff[0] if side in 'lr' else idxoff[1]
            idx, = np.where(panels==num)
            idx = slice(off + min(idx), off + max(idx) + 1)
            flush = kwargs[side + 'flush']
            if side=='r':
                subspec = gs[idx,-1]
                wspace, hspace = kwargs['rsep'], []
                wratios, hratios = kwargs['rwidth'], 1
                nrows, ncols = 1, len(wratios)
            elif side=='l':
                subspec = gs[idx,0]
                wspace, hspace = kwargs['lsep'], []
                wratios, hratios = kwargs['lwidth'], 1
                nrows, ncols = 1, len(wratios)
            elif side=='b':
                subspec = gs[-1,idx]
                wspace, hspace = [], kwargs['bsep']
                wratios, hratios = 1, kwargs['bwidth']
                nrows, ncols = len(hratios), 1
            # Make gridspec for containing the "stack" of panels
            ipaxs = []
            igs = gridspec.FlexibleGridSpecFromSubplotSpec(subplot_spec=subspec,
                    nrows=nrows, ncols=ncols,
                    wspace=wspace, hspace=hspace,
                    width_ratios=wratios, height_ratios=hratios,
                    )
            for i in range(max((nrows,ncols))):
                ipax = fig.add_subplot(igs[i], projection='panel', side=name, flush=flush)
                ipaxs += [ipax]
            paxs += [ipaxs]
        # Then *unfurl* lists of lists in consistent order
        if (side=='b' and order=='C') or (side in 'lr' and order=='F'):
            paxs = [*zip(*paxs)]
        paxs = axes_list([ax for ipaxs in paxs for ax in ipaxs], order=order, n=len(paxs[0]))
        setattr(fig, name + 'panel', paxs)

    # Set up axis sharing
    # NOTE: You must loop through twice, separately for x and y sharing. Fails
    # otherwise, don't know why.
    if sharex:
        for idx in range(naxs):
            igroup = np.where([idx in g for g in xgroups])[0]
            if igroup.size!=1:
                continue
            sharex_ax = axs[xgroups_base[igroup[0]]]
            if axs[idx] is sharex_ax:
                continue
            axs[idx]._sharex_setup(sharex_ax, sharex)
    if sharey:
        for idx in range(naxs):
            igroup = np.where([idx in g for g in ygroups])[0] # np.where works on lists
            if igroup.size!=1:
                continue
            sharey_ax = axs[ygroups_base[igroup[0]]]
            if axs[idx] is sharey_ax:
                continue
            axs[idx]._sharey_setup(sharey_ax, sharey)
    for ax in axs: # check axes don't belong to multiple groups, should be impossible unless my code is completely wrong...
        for name,groups in zip(('sharex', 'sharey'), (xgroups, ygroups)):
            if sum(ax in group for group in xgroups)>1:
                raise ValueError(f'Something went wrong; axis {idx:d} belongs to multiple {name} groups.')

    # Return results
    fig._main_axes = axs
    fig._ref_num = ref
    return fig, axes_list(axs)

