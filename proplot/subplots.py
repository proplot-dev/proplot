#!/usr/bin/env python3
"""
The starting point for creating custom ProPlot figures and axes.
The `subplots` command is all you'll need to directly use here;
it returns a figure instance and list of axes.

Note that
instead of separating features into their own functions (e.g.
a `generate_panel` function), we specify them right when the
figure is declared.

The reason for this approach is we want to build a
*static figure "scaffolding"* before plotting anything. This
way, ProPlot can hold the axes aspect ratios, panel widths, colorbar
widths, and inter-axes spaces **fixed** (note this was really hard
to do and is **really darn cool**!).

See `~FigureBase.smart_tight_layout` for details.
"""
import os
import re
import numpy as np
# Local modules, projection sand formatters and stuff
from .rcmod import rc
from .utils import _dot_dict, _default, _timer, _counter, units, ic
from . import gridspec, axes
import functools
import matplotlib.pyplot as plt
import matplotlib.figure as mfigure
import matplotlib.transforms as mtransforms
# For panel names
_aliases = {
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'tpanel': 'toppanel',
    'lpanel': 'leftpanel'
    }

#------------------------------------------------------------------------------#
# Miscellaneous helper functions
#------------------------------------------------------------------------------#
def close():
    """
    Alias for ``matplotlib.pyplot.close('all')``.
    Closes all figures stored in memory. Note this does not delete
    rendered figures in an iPython notebook.
    """
    plt.close('all') # easy peasy

def show():
    """
    Alias for ``matplotlib.pyplot.show()``. Note this command should
    be unnecessary if you are doing inline iPython notebook plotting
    and ran the `~proplot.notebook.nbsetup` command.
    """
    plt.show()

#------------------------------------------------------------------------------#
# Figure class
#------------------------------------------------------------------------------#
class FigureBase(mfigure.Figure):
    # Subclass adding some super cool features
    def __init__(self, figsize, gridspec=None, subplots_kw=None,
            tight=None, auto_adjust=True, pad=0.1,
            rcreset=True, silent=True, # print various draw statements?
            **kwargs):
        """
        Matplotlib figure with some pizzazz.

        Parameters
        ----------
        figsize : length-2 list of float or str
            Figure size (width, height). If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        gridspec : None or `~proplot.gridspec.FlexibleGridSpec`
            Required if `tight` is ``True``.
            The gridspec instance used to fill entire figure.
        subplots_kw : None or dict-like
            Required if `tight` is ``True``.
            Container of the keyword arguments used in the initial call
            to `subplots`.
        tight : bool, optional
            When figure is drawn, conform its size to a tight bounding
            box around figure contents, without messing
            up axes aspect ratios and internal spacing?
            See `~FigureBase.smart_tight_layout` for details.
        auto_adjust : bool, optional
            Alias for ``tight``.
        pad : float or str, optional
            Margin size around the tight bounding box. Ignored if
            `tight` is ``False``. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        rcreset : bool, optional
            Whether to reset all `~proplot.rcmod.rc` settings to their
            default values once the figure is drawn.
        silent : bool, optional
            Whether to silence messages related to drawing and printing
            the figure.

        Other parameters
        ----------------
        **kwargs : dict-like
            Passed to `matplotlib.figure.Figure` initializer.
        """
        # Initialize figure with some custom attributes.
        # Whether to reset rcParams wheenver a figure is drawn (e.g. after
        # ipython notebook finishes executing)
        self._rcreset  = rcreset
        self._smart_pad = units(pad)
        self._smart_tight = _default(tight, auto_adjust) # note name _tight already taken!
        self._smart_tight_init = True # is figure in its initial state?
        self._span_labels = [] # add axis instances to this, and label position will be updated
        self._silent = silent # whether to print message when we resize gridspec
        # Gridspec information
        self._gridspec = gridspec # gridspec encompassing drawing area
        self._subplots_kw = _dot_dict(subplots_kw) # extra special settings
        # Figure dimensions
        figsize = [units(size) for size in figsize]
        self.width  = figsize[0] # dimensions
        self.height = figsize[1]
        # Panels, initiate as empty
        self.leftpanel   = axes.EmptyPanel()
        self.bottompanel = axes.EmptyPanel()
        self.rightpanel  = axes.EmptyPanel()
        self.toppanel    = axes.EmptyPanel()
        # Proceed
        super().__init__(figsize=figsize, **kwargs) # python 3 only
        # Initialize suptitle, adds _suptitle attribute
        self.suptitle('')
        self._suptitle_init = True # whether vertical position has been changed yet
        self._suptitle_transform = None

    def __getattribute__(self, attr, *args):
        # Get attribute, but offer some convenient aliases
        attr = _aliases.get(attr, attr)
        return super().__getattribute__(attr, *args)

    def _rowlabels(self, labels, **kwargs):
        """
        Assign row labels.
        """
        axs = []
        for ax in self.axes:
            if isinstance(ax, axes.BaseAxes) and not isinstance(ax, axes.PanelAxes) and ax._col_span[0]==0:
                axs.append(ax)
        if isinstance(labels,str): # common during testing
            labels = [labels]*len(axs)
        if len(labels)!=len(axs):
            raise ValueError(f'Got {len(labels)} labels, but there are {len(axs)} rows.')
        axs = [ax for _,ax in sorted(zip([ax._row_span[0] for ax in axs],axs))]
        for ax,label in zip(axs,labels):
            if label and not ax.rowlabel.get_text():
                # Create a CompositeTransform that converts coordinates to
                # universal dots, then back to axes
                label_to_ax = ax.yaxis.label.get_transform() + ax.transAxes.inverted()
                x, _ = label_to_ax.transform(ax.yaxis.label.get_position())
                ax.rowlabel.set_visible(True)
                # Add text
                ax.rowlabel.update({'text':label,
                    'position':[x,0.5],
                    'ha':'right', 'va':'center', **kwargs})

    def _collabels(self, labels, **kwargs):
        """
        Assign column labels.
        """
        axs = []
        for ax in self.axes:
            if isinstance(ax, axes.BaseAxes) and not isinstance(ax, axes.PanelAxes) and ax._row_span[0]==0:
                axs.append(ax)
        if isinstance(labels,str):
            labels = [labels]*len(axs)
        if len(labels)!=len(axs):
            raise ValueError(f'Got {len(labels)} labels, but there are {len(axs)} columns.')
        axs = [ax for _,ax in sorted(zip([ax._col_span[0] for ax in axs],axs))]
        for ax,label in zip(axs,labels):
            if label and not ax.collabel.get_text():
                ax.collabel.update({'text':label, **kwargs})

    def _suptitle_setup(self, renderer=None, auto=True, **kwargs):
        """
        Intelligently determine super title position.

        Parameters
        ----------
        renderer : None or `~matplotlib.backend_bases.RendererBase`, optional
            The renderer.
        auto : bool, optional
            Whether this was called manually or not. In the latter case, we
            set a flag so that the vertical position will not be changed again.

        Notes
        -----
        Seems that `~matplotlib.figure.Figure.draw` is called *more than once*,
        and the title position are appropriately offset only during the
        *later* calls! So need to run this every time.
        """
        # Intelligently determine supertitle position:
        # Determine x by the underlying gridspec structure, where main axes lie.
        left = self._subplots_kw.left
        right = self._subplots_kw.right
        if self.leftpanel:
            left += (self._subplots_kw.lwidth + self._subplots_kw.lspace)
        if self.rightpanel:
            right += (self._subplots_kw.rwidth + self._subplots_kw.rspace)
        xpos = left/self.width + 0.5*(self.width - left - right)/self.width

        # if self._suptitle_init or kwargs.get('text', None):
        if False:
            # base = rc['axes.titlepad']/72 + self._gridspec.top*self.height
            # ypos = base/self.height
            ypos = self._suptitle.get_position()[1]
        else:
            # Figure out which title on the top-row axes will be offset the most
            # WARNING: Have to use private API to figure out whether axis has
            # tick labels or not! Seems to be no other way to do it.
            # See: https://matplotlib.org/_modules/matplotlib/axis.html#Axis.set_tick_params
            # TODO: Modify self.axes? Maybe add a axes_main and axes_panel
            # attribute or something? Because if have insets and other stuff
            # won't that fuck shit up?
            title_lev1, title_lev2, title_lev3 = None, None, None
            for ax in self.axes:
                if not isinstance(ax, axes.BaseAxes) or not ax._row_span[0]==0 or (isinstance(ax, axes.PanelAxes) and ax.panel_side=='bottom'):
                    continue
                title_lev1 = ax.title # always will be non-None
                if (ax.title.get_text() and not ax._title_inside) or ax.collabel.get_text():
                    title_lev2 = ax.title
                if ax.xaxis.get_ticks_position() == 'top':
                    test = 'label1On' not in ax.xaxis._major_tick_kw \
                        or ax.xaxis._major_tick_kw['label1On'] \
                        or ax.xaxis._major_tick_kw['label2On']
                    if test:
                        title_lev3 = ax.title

            # Hacky bugfixes
            # NOTE: Default linespacing is 1.2; it has no get, only a setter.
            # See: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
            if title_lev3: # upper axes present
                # If title and tick labels on top, offset the suptitle and get
                # matplotlib to adjust tight_subplot by prepending newlines to
                # the actual title.
                title = title_lev3
                text = title.get_text()
                title.set_text('\n\n' + text)
                line = 1.2 # looks best empirically
                add = '\n'
            elif title_lev2:
                # Just offset suptitle, and matplotlib will recognize the
                # suptitle during tight_subplot adjustment.
                line = 1.2
                title = title_lev2 # most common one
            else:
                # No title present. So, fill with spaces. Does nothing in most
                # cases, but if tick labels are on top, without this step
                # matplotlib tight subplots will not see the suptitle; now
                # suptitle will just occupy the empty title space.
                line = 0
                title = title_lev1
                if not title.axes._title_inside:
                    title.set_text(' ')
                add = ''

            # First idea: Create blended transform, end with newline
            # ypos = title.axes._title_pos_init[1]
            # kwargs['text'] = kwargs.get('text', self._suptitle.get_text()) + add
            # kwargs['transform'] = mtransforms.blended_transform_factory(self.transFigure, title.get_transform())

            # New idea: Get the transformed position
            if self._suptitle_transform:
                transform = self._suptitle_transform
            else:
                transform = title.axes._title_pos_transform + self.transFigure.inverted()
            ypos = transform.transform(title.axes._title_pos_init)[1]
            line = line*(title.get_size()/72)/self.height
            ypos = ypos + line
            kwargs['transform'] = self.transFigure

        # Update settings
        self._suptitle.update({'position':(xpos, ypos),
            'ha':'center', 'va':'bottom', **kwargs})
        if auto:
            self._suptitle_init = False

    def draw(self, renderer, *args, **kwargs):
        """
        Fix the "super title" position and automatically adjust the
        main gridspec, then call the parent `~matplotlib.figure.Figure.draw`
        method.
        """
        # Special: Figure out if other titles are present, and if not
        # bring suptitle close to center
        self._suptitle_setup(renderer, auto=True) # just applies the spacing
        self._auto_smart_tight_layout(renderer)
        out = super().draw(renderer, *args, **kwargs)

        # If rc settings have been changed, reset them after drawing is done
        # Usually means we have finished executing a notebook cell
        if not rc._init and self._rcreset:
            if not self._silent:
                print('Resetting rcparams.')
            rc.reset()
        return out

    def panel_factory(self, subspec, whichpanels=None,
            hspace=None, wspace=None,
            hwidth=None, wwidth=None,
            sharex=None, sharey=None, # external sharing
            sharex_level=3, sharey_level=3,
            sharex_panels=True, sharey_panels=True, # by default share main x/y axes with panel x/y axes
            **kwargs):
        """
        Create "inner" panels, i.e. panels along subplot edges.

        Parameters
        ----------
        subspec : `~matplotlib.gridspec.SubplotSpec`
            The `~matplotlib.gridspec.SubplotSpec` instance onto which
            the main subplot and its panels are drawn.
        whichpanels : None or str, optional
            Whether to draw panels on the left, right, bottom, or top
            sides. Should be a string containing any of the characters
            ``'l'``, ``'r'``, ``'b'``, or ``'t'`` in any order.
            Default is ``'r'``.
        wwidth, hwidth : None, float, or str, optional
            Width of vertical (`wwidth`) and horizontal (`hwidth`)
            panels, respectively. If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        wspace, hspace : None, float, or str, optional
            Empty space between the main subplot and the panel.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`.
        sharex_panels, sharey_panels : bool, optional
            Whether to share the *x* axes labels for vertically oriented
            (left/right) panels stacked on top of each other, and the
            *y* axes labels for horizontally oriented (bottom/top) panels
            next to each other. Defaults to ``True``.
        sharex, sharey : None or `~proplot.axes.BaseAxes`, optional
            The axes for "sharing" labels, tick labels, and limits.
            See `subplots` for details.
        sharex_level, sharey_level : {3, 2, 1, 0}, optional
            The axis sharing level.
            See `subplots` for details.

        Todo
        ----
        Make settings specific to left, right, top, bottom sides!
        """
        # Helper function for creating paneled axes.
        width, height = self.width, self.height
        translate = {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}
        whichpanels = translate.get(whichpanels, whichpanels)
        whichpanels = whichpanels or 'r'
        hspace = units(_default(hspace, 0.13)) # teeny tiny space
        wspace = units(_default(wspace, 0.13))
        hwidth = units(_default(hwidth, 0.45)) # default is panels for plotting stuff, not colorbars
        wwidth = units(_default(wwidth, 0.45))
        if any(s.lower() not in 'lrbt' for s in whichpanels):
            raise ValueError(f'Whichpanels argument can contain characters l (left), r (right), b (bottom), or t (top), instead got "{whichpanels}".')

        # Determine rows/columns and indices
        nrows = 1 + sum(1 for i in whichpanels if i in 'bt')
        ncols = 1 + sum(1 for i in whichpanels if i in 'lr')
        sides_lr = [l for l in ['l',None,'r'] if not l or l in whichpanels]
        sides_tb = [l for l in ['t',None,'b'] if not l or l in whichpanels]
        # Detect empty positions and main axes position
        main_pos  = (int('t' in whichpanels), int('l' in whichpanels))
        corners   = {'tl':(0,0),             'tr':(0,main_pos[1]+1),
                     'bl':(main_pos[0]+1,0), 'br':(main_pos[0]+1,main_pos[1]+1)}
        empty_pos = [position for corner,position in corners.items() if
                     corner[0] in whichpanels and corner[1] in whichpanels]

        # Fix wspace/hspace in inches, using the Bbox from get_postition
        # on the subspec object to determine physical width of axes to be created
        # * Consider writing some convenience funcs to automate this unit conversion
        bbox = subspec.get_position(self) # valid since axes not drawn yet
        if hspace is not None:
            hspace = np.atleast_1d(hspace)
            if hspace.size==1:
                hspace = np.repeat(hspace, (nrows-1,))
            boxheight = np.diff(bbox.intervaly)[0]*height
            height = boxheight - hspace.sum()
            hspace = hspace/(height/nrows)
        if wspace is not None:
            wspace = np.atleast_1d(wspace)
            if wspace.size==1:
                wspace = np.repeat(wspace, (ncols-1,))
            boxwidth = np.diff(bbox.intervalx)[0]*width
            width = boxwidth - wspace.sum()
            wspace = wspace/(width/ncols)

        # Figure out hratios/wratios
        # Will enforce (main_width + panel_width)/total_width = 1
        wwidth_ratios = [width - wwidth*(ncols-1)]*ncols
        if wwidth_ratios[0]<0:
            raise ValueError(f'Panel wwidth {wwidth} is too large. Must be less than {width/(nrows-1):.3f}.')
        for i in range(ncols):
            if i!=main_pos[1]: # this is a panel entry
                wwidth_ratios[i] = wwidth
        hwidth_ratios = [height - hwidth*(nrows-1)]*nrows
        if hwidth_ratios[0]<0:
            raise ValueError(f'Panel hwidth {hwidth} is too large. Must be less than {height/(ncols-1):.3f}.')
        for i in range(nrows):
            if i!=main_pos[0]: # this is a panel entry
                hwidth_ratios[i] = hwidth

        # Create subplotspec and draw the axes
        # Will create axes in order of rows/columns so that the "base" axes
        # are always built before the axes to be "shared" with them
        panels = []
        gs = gridspec.FlexibleGridSpecFromSubplotSpec(
                nrows         = nrows,
                ncols         = ncols,
                subplot_spec  = subspec,
                wspace        = wspace,
                hspace        = hspace,
                width_ratios  = wwidth_ratios,
                height_ratios = hwidth_ratios,
                )
        # Draw main axes
        ax = self.add_subplot(gs[main_pos[0], main_pos[1]], **kwargs)
        axmain = ax
        # Draw axes
        panels = {}
        kwpanels = {**kwargs, 'projection':'panel'} # override projection
        kwpanels.pop('number', None) # don't want numbering on panels
        translate = {'b':'bottom', 't':'top', 'l':'left', 'r':'right'} # inverse
        for r,side_tb in enumerate(sides_tb): # iterate top-bottom
            for c,side_lr in enumerate(sides_lr): # iterate left-right
                if (r,c) in empty_pos or (r,c)==main_pos:
                    continue
                side = translate.get(side_tb or side_lr, None)
                ax = self.add_subplot(gs[r,c], panel_side=side, panel_parent=axmain, **kwpanels)
                panels[side] = ax

        # Finally add as attributes, and set up axes sharing
        axmain.bottompanel = panels.get('bottom', axes.EmptyPanel())
        axmain.toppanel    = panels.get('top',    axes.EmptyPanel())
        axmain.leftpanel   = panels.get('left',   axes.EmptyPanel())
        axmain.rightpanel  = panels.get('right',  axes.EmptyPanel())
        if sharex_panels:
            axmain._sharex_panels()
        if sharey_panels:
            axmain._sharey_panels()
        axmain._sharex_setup(sharex, sharex_level)
        axmain._sharey_setup(sharey, sharey_level)
        return axmain

    def smart_tight_layout(self, renderer=None):
        """
        Conform figure edges to tight bounding box around content, without
        screwing up subplot aspect ratios, empty spaces, panel sizes, and such.

        This is called automatically whenever `~FigureBase.draw` or
        `~FigureBase.savefig` are called,
        unless the user set `auto_adjust` to ``False`` in the call to
        `subplots`.

        Parameters
        ----------
        renderer : None or `~matplotlib.backend_bases.RendererBase`, optional
            The backend renderer.
        """
        # Get bounding box that encompasses *all artists*, compare to bounding
        # box used for saving *figure*
        pad = self._smart_pad
        if self._subplots_kw is None or self._gridspec is None:
            raise ValueError("Initialize figure with 'subplots_kw' and 'gridspec' to draw tight grid.")
        obbox = self.bbox_inches # original bbox
        if not renderer: # cannot use the below on figure save! figure becomes a special FigurePDF class or something
            renderer = self.canvas.get_renderer()
        bbox = self.get_tightbbox(renderer)
        ox, oy, x, y = obbox.intervalx, obbox.intervaly, bbox.intervalx, bbox.intervaly
        x1, y1, x2, y2 = x[0], y[0], ox[1]-x[1], oy[1]-y[1] # deltas

        # Apply new settings
        lname = 'lspace' if self.leftpanel else 'left'
        rname = 'rspace' if self.rightpanel else 'right'
        bname = 'bspace' if self.bottompanel else 'bottom'
        tname = 'top'
        subplots_kw = self._subplots_kw
        left   = getattr(subplots_kw, lname) - x1 + pad
        right  = getattr(subplots_kw, rname) - x2 + pad
        bottom = getattr(subplots_kw, bname) - y1 + pad
        top    = getattr(subplots_kw, tname) - y2 + pad
        subplots_kw.update({lname:left, rname:right, bname:bottom, tname:top})
        figsize, *_, gridspec_kw = _subplots_kwargs(**subplots_kw)
        self._smart_tight_init = False
        self._gridspec.update(**gridspec_kw)
        self.set_size_inches(figsize)

        # Fix any spanning labels that we've added to _span_labels
        # These need figure-relative height coordinates
        for axis in self._span_labels:
            axis.axes._share_span_label(axis)

    def _auto_smart_tight_layout(self, renderer=None):
        # Only proceed if this wasn't already done, or user did not want
        # the figure to have tight boundaries.
        if not self._smart_tight_init or not self._smart_tight:
            return
        # Cartopy sucks at labels! Bounding box identified will be wrong.
        # 1) If you used set_bounds to zoom into part of a cartopy projection,
        # this can erroneously identify invisible edges of map as being part of boundary
        # 2) If you have gridliner text labels, matplotlib won't detect them.
        # Therefore, bail if we find any cartopy axes.
        # TODO: Fix this?
        if any(isinstance(ax, axes.CartopyAxes) for ax in self.axes):
            return
        # Proceed
        if not self._silent:
            print('Adjusting gridspec.')
        self.smart_tight_layout(renderer)

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
        # Notes:
        # * Gridspec object must be updated before figure is printed to
        #     screen in interactive environment; will fail to update after that.
        #     Seems to be glitch, should open thread on GitHub.
        # * To color axes patches, you may have to explicitly pass the
        #     transparent=False kwarg.
        #     Some kwarg translations, to pass to savefig
        if 'alpha' in kwargs:
            kwargs['transparent'] = not bool(kwargs.pop('alpha')) # 1 is non-transparent
        if 'color' in kwargs:
            kwargs['facecolor'] = kwargs.pop('color') # the color
            kwargs['transparent'] = False
        # Finally, save
        self._suptitle_setup(auto=True) # just applies the spacing
        self._auto_smart_tight_layout()
        if not self._silent:
            print(f'Saving to "{filename}".')
        return super().savefig(os.path.expanduser(filename), **kwargs) # specify DPI for embedded raster objects

    def save(self, *args, **kwargs):
        """
        Alias for `~FigureBase.savefig`.
        """
        # Alias for save.
        return self.savefig(*args, **kwargs)

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
class axes_list(list):
    """
    Magical class that iterates through items and calls respective
    method (or retrieves respective attribute) on each one. For example,
    ``axs.format(color='r')`` colors all axes spines and tick marks red.

    When using `subplots`, the list of axes returned is an instance of
    `axes_list`.
    """
    def __repr__(self):
        # Make clear that this is no ordinary list
        return 'axes_list(' + super().__repr__() + ')'

    def __getitem__(self, key):
        # Return an axes_list version of the slice, or just the axes
        axs = list.__getitem__(self, key)
        if isinstance(key,slice): # i.e. returns a list
            axs = axes_list(axs)
        return axs

    def __getattr__(self, attr):
        # Stealthily return dummy function that actually iterates
        # through each attribute here
        values = [getattr(ax, attr, None) for ax in self]
        if None in values:
            raise AttributeError(f"'{type(self[0])}' object has no method '{attr}'.")
        elif all(callable(value) for value in values):
            @functools.wraps(values[0])
            def iterator(*args, **kwargs):
                ret = []
                for ax in self:
                    res = getattr(ax, attr)(*args, **kwargs)
                    if res is not None:
                        ret += [res]
                return None if not ret else ret[0] if len(ret)==1 else ret
            return iterator
        elif all(not callable(value) for value in values):
            return values[0] if len(values)==1 else values # just return the attribute list
        else:
            raise AttributeError('Mixed methods found.')

# Function for processing input and generating necessary keyword args
def _subplots_kwargs(nrows, ncols, rowmajor=True,
    aspect=1,    figsize=None, # for controlling aspect ratio, default is control for width
    # General
    width=None,  height=None, axwidth=None, axheight=None, journal=None,
    hspace=None, wspace=None, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
    left=None,   bottom=None, right=None,   top=None,     # spaces around edge of main plotting area, in inches
    bwidth=None, bspace=None, rwidth=None, rspace=None, lwidth=None, lspace=None, # default to no space between panels
    # Panels
    bottompanel=False, bottompanels=False,
    rightpanel=False,  rightpanels=False,
    leftpanel=False,   leftpanels=False,
    bottomcolorbar=False, bottomcolorbars=False, bottomlegend=False, bottomlegends=False, # convenient aliases that change default features
    rightcolorbar=False,  rightcolorbars=False,  rightlegend=False,  rightlegends=False,
    leftcolorbar=False,   leftcolorbars=False,   leftlegend=False,   leftlegends=False,
    # Panel aliases
    bpanel=None,    bpanels=None,
    rpanel=None,    rpanels=None,
    lpanel=None,    lpanels=None,
    bcolorbar=None, bcolorbars=None, blegend=None, blegends=None,
    rcolorbar=None, rcolorbars=None, rlegend=None, rlegends=None,
    lcolorbar=None, lcolorbars=None, llegend=None, llegends=None,
    ):
    """
    Handle complex keyword args and aliases thereof.
    """
    # Aliases
    # Because why the fuck not?
    bottompanel  = bpanel or bottompanel
    bottompanels = bpanels or bottompanels
    rightpanel   = rpanel or rightpanel
    rightpanels  = rpanels or rightpanels
    leftpanel    = lpanel or leftpanel
    leftpanels   = lpanels or leftpanels
    bottomlegend  = blegend or bottomlegend
    bottomlegends = blegends or bottomlegends # convenient aliases that change default features
    rightlegend   = rlegend or rightlegend
    rightlegends  = rlegends or rightlegends
    leftlegend    = llegend or leftlegend
    leftlegends   = llegends or leftlegends
    bottomcolorbar  = bcolorbar or bottomcolorbar
    bottomcolorbars = bcolorbars or bottomcolorbars
    rightcolorbar   = rcolorbar or rightcolorbar
    rightcolorbars  = rcolorbars or rightcolorbars
    leftcolorbar    = lcolorbar or leftcolorbar
    leftcolorbars   = lcolorbars or leftcolorbars

    # Handle the convenience feature for specifying the desired width/spacing
    # for panels as that suitable for a colorbar or legend
    # NOTE: Ugly but this is mostly boilerplate, shouln't change much
    def _panelprops(panel, panels, colorbar, colorbars, legend, legends, width, space):
        if colorbar or colorbars:
            width = _default(width, rc['gridspec.cbar'])
            space = _default(space, rc['gridspec.xlab'])
            panel, panels = colorbar, colorbars
        elif legend or legends:
            width = _default(width, rc['gridspec.legend'])
            space = _default(space, 0)
            panel, panels = legend, legends
        return panel, panels, width, space
    rightpanel, rightpanels, rwidth, rspace, = _panelprops(
        rightpanel, rightpanels, rightcolorbar, rightcolorbars,
        rightlegend, rightlegends, rwidth, rspace)
    leftpanel, leftpanels, lwidth, lspace = _panelprops(
        leftpanel, leftpanels, leftcolorbar, leftcolorbars,
        leftlegend, leftlegends, lwidth, lspace)
    bottompanel, bottompanels, bwidth, bspace = _panelprops(
        bottompanel, bottompanels, bottomcolorbar, bottomcolorbars,
        bottomlegend, bottomlegends, bwidth, bspace)

    # Handle the convenience feature for generating one panel per row/column
    # and one single panel for all rows/columns
    def _parse(panel, panels, nmax):
        if panel: # one spanning panel
            panels = [1]*nmax
        elif panels not in (None,False): # can't test truthiness, want user to be allowed to pass numpy vector!
            try:
                panels = list(panels)
            except TypeError:
                panels = [*range(nmax)] # pass True to make panel for each column
        return panels
    bottompanels = _parse(bottompanel, bottompanels, ncols)
    rightpanels  = _parse(rightpanel,  rightpanels,  nrows)
    leftpanels   = _parse(leftpanel,   leftpanels,   nrows)

    # Apply the general defaults
    # Need to do this after number of rows/columns figured out
    wratios = np.atleast_1d(_default(wratios, 1))
    hratios = np.atleast_1d(_default(hratios, 1))
    hspace = np.atleast_1d(_default(hspace, rc['gridspec.title']))
    wspace = np.atleast_1d(_default(wspace, rc['gridspec.inner']))
    if len(wratios)==1:
        wratios = np.repeat(wratios, (ncols,))
    if len(hratios)==1:
        hratios = np.repeat(hratios, (nrows,))
    if len(wspace)==1:
        wspace = np.repeat(wspace, (ncols-1,))
    if len(hspace)==1:
        hspace = np.repeat(hspace, (nrows-1,))
    left   = units(_default(left,   rc['gridspec.ylab']))
    bottom = units(_default(bottom, rc['gridspec.xlab']))
    right  = units(_default(right,  rc['gridspec.nolab']))
    top    = units(_default(top,    rc['gridspec.title']))
    bwidth = units(_default(bwidth, rc['gridspec.cbar']))
    rwidth = units(_default(rwidth, rc['gridspec.cbar']))
    lwidth = units(_default(lwidth, rc['gridspec.cbar']))
    bspace = units(_default(bspace, rc['gridspec.xlab']))
    rspace = units(_default(rspace, rc['gridspec.ylab']))
    lspace = units(_default(lspace, rc['gridspec.ylab']))

    # Determine figure size
    if journal:
        if width or height or axwidth or axheight or figsize:
            raise ValueError('Argument conflict: Specify only a journal size, or the figure dimensions, not both.')
        width, height = journal_size(journal) # if user passed width=<string>, will use that journal size
    if not figsize:
        figsize = (width, height)
    width, height = figsize
    width  = units(width, error=False)
    height = units(height, error=False)

    # If width and height are not fixed, determine necessary width/height to
    # preserve the aspect ratio of specified plot
    auto_both = (width is None and height is None)
    auto_width  = (width is None and height is not None)
    auto_height = (height is None and width is not None)
    auto_neither = (width is not None and height is not None)
    bpanel_space = bwidth + bspace if bottompanels else 0
    rpanel_space = rwidth + rspace if rightpanels else 0
    lpanel_space = lwidth + lspace if leftpanels else 0
    try:
        aspect = aspect[0]/aspect[1]
    except (IndexError,TypeError):
        pass # do nothing
    aspect_fixed = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    aspect_fixed = aspect*(hratios[0]/np.mean(hratios))
    # Determine average axes widths/heights
    # Default behavior: axes average 2.0 inches wide
    if auto_width or auto_neither:
        axheight_ave = (height - top - bottom - sum(hspace) - bpanel_space)/nrows
    if auto_height or auto_neither:
        axwidth_ave = (width - left - right - sum(wspace) - rpanel_space - lpanel_space)/ncols
    if auto_both: # get stuff directly from axes
        if axwidth is None and axheight is None:
            axwidth = 2.0
        if axheight is not None:
            height = axheight*nrows + top + bottom + sum(hspace) + bpanel_space
            auto_width = True
            axheight_ave = axheight
        if axwidth is not None:
            width = axwidth*ncols + left + right + sum(wspace) + rpanel_space + lpanel_space
            auto_height = True
            axwidth_ave = axwidth
        if axwidth is not None and axheight is not None:
            auto_width = auto_height = False
        figsize = (width, height) # again
    # Fix height and top-left axes aspect ratio
    if auto_width:
        axwidth_ave = axheight_ave*aspect_fixed
        width       = axwidth_ave*ncols + left + right + sum(wspace) + rpanel_space + lpanel_space
    # Fix width and top-left axes aspect ratio
    if auto_height:
        axheight_ave = axwidth_ave/aspect_fixed
        height       = axheight_ave*nrows + top + bottom + sum(hspace) + bpanel_space
    # Check
    if axwidth_ave<0:
        raise ValueError(f"Not enough room for axes (would have width {axwidth_ave}). Increase width, or reduce spacings 'left', 'right', or 'wspace'.")
    if axheight_ave<0:
        raise ValueError(f"Not enough room for axes (would have height {axheight_ave}). Increase height, or reduce spacings 'top', 'bottom', or 'hspace'.")

    # Necessary arguments to reconstruct this grid
    # Can follow some of the pre-processing
    subplots_kw = _dot_dict(nrows=nrows, ncols=ncols,
        figsize=figsize, aspect=aspect,
        hspace=hspace,   wspace=wspace,
        hratios=hratios, wratios=wratios,
        bottompanels=bottompanels, leftpanels=leftpanels, rightpanels=rightpanels,
        left=left,     bottom=bottom, right=right,   top=top,
        bwidth=bwidth, bspace=bspace, rwidth=rwidth, rspace=rspace, lwidth=lwidth, lspace=lspace,
        )

    # Make sure the 'ratios' and 'spaces' are in physical units (we cast the
    # former to physical units), easier then to add stuff as below
    wspace = wspace.tolist()
    hspace = hspace.tolist()
    wratios = (ncols*axwidth_ave*(wratios/sum(wratios))).tolist()
    hratios = (nrows*axheight_ave*(hratios/sum(hratios))).tolist()

    # Now add the outer panel considerations (idea is we have panels whose
    # widths/heights are *in inches*, and only allow the main subplots and
    # figure widhts/heights to warp to preserve aspect ratio)
    nrows += int(bool(bottompanels))
    ncols += int(bool(rightpanels)) + int(bool(leftpanels))
    if bottompanels: # the 'bottom' space actually goes between subplots and panel
        hratios = hratios + [bwidth] # easy
        hspace  = hspace + [bottom]
        bottom  = bspace
    if leftpanels:
        wratios = [lwidth] + wratios
        wspace  = [left] + wspace
        left    = lspace
    if rightpanels:
        wratios = wratios + [rwidth]
        wspace  = wspace + [right]
        right   = rspace
    # Scale stuff that gridspec needs to be scaled
    # Scale the boundaries for gridspec
    # NOTE: We *no longer* scale wspace/hspace because we expect it to
    # be in same scale as axes ratios, much easier that way and no drawback really
    bottom = bottom/height
    left   = left/width
    top    = 1-top/height
    right  = 1-right/width

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    offset = (0, 1 if leftpanels else 0)
    figsize = (width, height)
    gridspec_kw = dict(
            nrows         = nrows,
            ncols         = ncols,
            left          = left,
            bottom        = bottom,
            right         = right, # unique spacing considerations
            top           = top, # so far no panels allowed here
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = wratios,
            height_ratios = hratios,
            ) # set wspace/hspace to match the top/bottom spaces
    return figsize, offset, subplots_kw, gridspec_kw

def subplots(array=None, ncols=1, nrows=1,
        order='C', # allow calling with subplots(array)
        emptycols=[], emptyrows=[], # obsolete?
        tight=None, auto_adjust=True, pad=0.1,
        rcreset=True, silent=True, # arguments for figure instantiation
        span=None, # bulk apply to x/y axes
        share=None, # bulk apply to x/y axes
        spanx=1,  spany=1,  # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        sharex=3, sharey=3, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        innerpanels={}, innercolorbars={}, innerpanels_kw={}, innercolorbars_kw={},
        basemap=False, proj={}, projection={}, proj_kw={}, projection_kw={},
        **kwargs):
    # See `proplot.axes.XYAxes.format`, `proplot.axes.CartopyAxes.format`, and
    # `proplot.axes.BasemapAxes.format` for formatting your figure.
    """
    Analagous to `matplotlib.pyplot.subplots`. Create a figure with a single
    axes or arbitrary grids of axes (any of which can be map projections),
    and optional "panels" along axes or figure edges.

    The parameters are sorted into the following rough sections: subplot grid
    specifications, figure and subplot sizes, axis sharing,
    inner panels, outer panels, and map projections.

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
    emptyrows, emptycols : list of int, optional
        Row, column numbers (starting from 1) that you want empty. Generally
        this is used with `ncols` and `nrows`; with `array`, you can
        just use zeros.

    width, height : float or str, optional
        The figure width and height. If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units` -- for example,
        `width="10cm"` creates a 10cm wide figure.
    figsize : length-2 tuple, optional
        Tuple specifying the figure `(width, height)`.
    journal : None or str, optional
        Conform figure width (and height, if specified) to academic journal
        standards. See `~proplot.gridspec.journal_size` for details.
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
    hspace, wspace, hratios, wratios, left, right, top, bottom : optional
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

    innerpanels : str or dict-like, optional
        Specify which axes should have inner "panels".

        Panels are stored on the ``leftpanel``, ``rightpanel``,
        ``bottompanel`` and ``toppanel`` attributes on the axes instance.
        You can also use the aliases ``lpanel``, ``rpanel``, ``bpanel``,
        or ``tpanel``.

        The argument is interpreted as follows:

            * If str, panels are drawn on the same side for all subplots.
              String should contain any of the characters ``'l'`` (left panel),
              ``'r'`` (right panel), ``'t'`` (top panel), or ``'b'`` (bottom panel).
              For example, ``'rt'`` indicates a right and top panel.
            * If dict-like, panels can be drawn on different sides for
              different subplots. For example, `innerpanels={1:'r', (2,3):'l'}`
              indicates that we want to draw a panel on the right side of
              subplot number 1, on the left side of subplots 2 and 3.

    innercolorbars : str or dict-like, optional
        Identical to ``innerpanels``, except the default panel size is more
        appropriate for a colorbar. The panel can then be **filled** with
        a colorbar with, for example, ``ax.leftpanel.colorbar(...)``.
    innerpanels_kw : dict-like, optional
        Keyword args passed to `~FigureBase.panel_factory`. Controls
        the width/height and spacing of panels.

        Can be dict of properties (applies globally), or **dict of dicts** of
        properties (applies to specific properties, as with `innerpanels`).

        For example, consider a figure with 2 columns and 1 row.
        With ``{'wwidth':1}``, all left/right panels 1 inch wide, while
        with ``{1:{'wwidth':1}, 2:{'wwidth':0.5}}``, the left subplot
        panel is 1 inch wide and the right subplot panel is 0.5 inches wide.

        See `~FigureBase.panel_factory` for keyword arg options.
    innercolorbars_kw
        Alias for `innerpanels_kw`.

    bottompanel, rightpanel, leftpanel : bool, optional
        Whether to draw an "outer" panel surrounding the entire grid of
        subplots, on the bottom, right, or left sides of the figure.
        This is great for creating global colorbars and legends.
        Default is ``False``.

        Panels are stored on the ``leftpanel``, ``rightpanel``,
        and ``bottompanel`` attributes on the figure instance.
        You can also use the aliases ``lpanel``, ``rpanel``, and ``bpanel``.
    bottompanels, rightpanels, leftpanels : bool or list of int, optional
        As with the `panel` keyword args, but this allows you to
        specify **multiple** panels along each figure edge.
        Default is ``False``.

        The arguments are interpreted as follows:

            * If bool, separate panels **for each column/row of subplots**
              are drawn.
            * If list of int, can specify panels that span **contiguous
              columns/rows**. Usage is similar to the ``array`` argument.

              For example, for a figure with 3 columns, ``bottompanels=[1, 2, 2]``
              draws a panel on the bottom of the first column and spanning the
              bottom of the right 2 columns, and ``bottompanels=[0, 2, 2]``
              only draws a panel underneath the right 2 columns (as with
              `array`, the ``0`` indicates an empty space).

    bottomcolorbar, bottomcolorbars, rightcolorbar,  rightcolorbars, leftcolorbar, leftcolorbars
        Identical to the corresponding ``panel`` and ``panels`` keyword
        args, except the default panel size is more appropriate for a colorbar.

        The panel can then be **filled** with a colorbar with, for example, 
        ``fig.leftpanel.colorbar(...)``.
    bottomlegend, bottomlegends, rightlegend,  rightlegends, leftlegend, leftlegends
        Identical to the corresponding ``panel`` and ``panels`` keyword
        args, except the default panel size is more appropriate for a legend.

        The panel can then be **filled** with a legend with, for example, 
        ``fig.leftpanel.legend(...)``.
    bpanel, bpanels, bcolorbar, bcolorbars, blegend, blegends, rpanel, rpanels, rcolorbar, rcolorbars, rlegend, rlegends, lpanel, lpanels, lcolorbar, lcolorbars, llegend, llegends
        Aliases for equivalent args with ``bottom``, ``right``, and ``left``.
        It's a bit faster, so why not?
    lwidth, rwidth, bwidth, twidth : float, optional
        Width of left, right, bottom, and top panels. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.
    lspace, rspace, bspace, tspace : float, optional
        Space between the "inner" edges of the left, right, bottom, and top
        panels and the edge of the main subplot grid. If float, units are
        inches. If string, units are interpreted by `~proplot.utils.units`.

    projection : str or dict-like, optional
        The map projection name.

        If string, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `innerpanels`.

        For example, consider a figure with 4 columns and 1 row.
        With ``projection={1:'mercator', (2,3):'hammer'}``,
        the leftmost subplot is a Mercator projection, the middle 2 are Hammer
        projections, and the rightmost is a normal x/y axes.
    projection_kw : dict-like, optional
        Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
        cartopy `~cartopy.crs.Projection` class on instantiation.

        As with `innerpanels_kw`,
        can be dict of properties (applies globally), or **dict of dicts** of
        properties (applies to specific properties, as with `innerpanels`).
    proj, proj_kw
        Aliases for `projection`, `projection_kw`.
    basemap : bool or dict-like, optional
        Whether to use `~mpl_toolkits.basemap.Basemap` or
        `~cartopy.crs.Projection` for map projections. Defaults to ``False``.

        If boolean, applies to all subplots. If dictionary, can be
        specific to each subplot, as with `innerpanels`.

    Returns
    -------
    f : `FigureBase`
        The figure instance.
    axs : `axes_list`
        A special list of axes instances. See `axes_list`.

    Other parameters
    ----------------
    auto_adjust, tight, pad, rcreset, silent
        Passed to `FigureBase`.
    **kwargs
        Passed to `~proplot.gridspec.FlexibleGridSpec`.

    See also
    --------
    `~proplot.axes.BaseAxes`, `~proplot.axes.BaseAxes.format`

    Notes
    -----
    * Matplotlib `~matplotlib.axes.Axes.set_aspect` option seems to behave
      strangely for some plots; for this reason we override the ``fix_aspect``
      keyword arg provided by `~mpl_toolkits.basemap.Basemap` and just draw
      the figure with appropriate aspect ratio to begin with. Otherwise, we
      get weird differently-shaped subplots that seem to make no sense.

    * All shared axes will generally end up with the same axis
      limits, scaling, major locators, and minor locators. The
      ``sharex`` and ``sharey`` detection algorithm really is just to get
      instructions to make the tick labels and axis labels **invisible**
      for certain axes.

    Todo
    ----
    * Add options for e.g. `bpanel` keyword args.
    * Fix axes aspect ratio stuff when width/height ratios are not one!
    * Generalize axes sharing for right y-axes and top x-axes. Enable a secondary
      axes sharing mode where we *disable ticklabels and labels*, but *do not
      use the builtin sharex/sharey API*, suitable for complex map projections.
    * For spanning axes labels, right now only detect **x labels on bottom**
      and **ylabels on top**. Generalize for all subplot edges.
    """
    # Check
    sharex = _default(share, sharex)
    sharey = _default(share, sharey)
    spanx  = _default(span, spanx)
    spany  = _default(span, spany)
    if int(sharex) not in range(4) or int(sharey) not in range(4):
        raise ValueError('Axis sharing options sharex/sharey can be 0 (no sharing), 1 (sharing, but keep all tick labels), and 2 (sharing, but only keep one set of tick labels).')
    # Helper functions
    translate = lambda p: {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}.get(p, p)
    auto_adjust = _default(tight, auto_adjust)
    def axes_dict(value, kw=False):
        # First build up dictionary
        # Accepts:
        # 1) 'string' or {1:'string1', (2,3):'string2'}
        if not kw:
            if not isinstance(value, dict):
                value = {range(1,num_axes+1): value}
        # 2) {'prop':value} or {1:{'prop':value1}, (2,3):{'prop':value2}}
        else:
            nested = [isinstance(value,dict) for value in value.values()]
            if not any(nested): # any([]) == False
                value = {range(1,num_axes+1): value.copy()}
            elif not all(nested):
                raise ValueError('Wut.')
        # Then unfurl wherever keys contain multiple axes numbers
        kw_out = {}
        for nums,item in value.items():
            nums = np.atleast_1d(nums)
            for num in nums.flat:
                if not kw:
                    kw_out[num-1] = item
                else:
                    kw_out[num-1] = item.copy()
        # Verify numbers
        if {*range(num_axes)} != {*kw_out.keys()}:
            raise ValueError(f'Have {num_axes} axes, but {value} has properties for axes {", ".join(str(i+1) for i in sorted(kw_out.keys()))}.')
        return kw_out

    # Array setup
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
    num_axes = len(nums)
    if tuple(nums.flat) != tuple(range(1,num_axes+1)):
        raise ValueError('Axes numbers must span integers 1 to num_axes (i.e. cannot skip over numbers).')
    nrows = array.shape[0]
    ncols = array.shape[1]

    # Get basemap.Basemap or cartopy.CRS instances for map, and override aspec tratio
    # NOTE: Previously went to some pains (mainly for basemap, something in the
    # initialization deals with this) to only draw one projection. This is hard
    # to generalize when want different projections/kwargs, so abandon
    basemap = axes_dict(basemap, False) # package used for projection
    proj = axes_dict(projection or proj or 'xy', False) # name of projection; by default use base.XYAxes
    proj_kw = axes_dict(projection_kw or proj_kw, True) # stores cartopy/basemap arguments
    axes_kw = {num:{} for num in range(num_axes)} # stores add_subplot arguments
    for num,name in proj.items():
        # Builtin matplotlib polar axes, just use my overridden version
        if name=='polar':
            axes_kw[num]['projection'] = 'newpolar'
            if num==1:
                kwargs.update(aspect=1)
        # The default, my XYAxes projection
        elif name=='xy':
            axes_kw[num]['projection'] = 'xy'
        # Custom Basemap and Cartopy axes
        elif name:
            package = 'basemap' if basemap[num] else 'cartopy'
            instance, aspect = axes.map_projection_factory(name, basemap=basemap[num], **proj_kw[num])
            axes_kw[num].update({'projection':package, 'map_projection':instance})
            if num==0:
                # print(f'Forcing aspect ratio: {aspect:.3g}')
                kwargs.update(aspect=aspect)
        else:
            raise ValueError('All projection names should be declared. Wut.')

    # Create dictionary of panel toggles and settings
    # Input can be string e.g. 'rl' or dictionary e.g. {(1,2,3):'r', 4:'l'}
    # NOTE: Internally we convert array references to 0-base here
    # Add kwargs and the 'which' arguments
    # Optionally change the default panel widths for 'colorbar' panels
    if not isinstance(innercolorbars, (dict, str)):
        raise ValueError('Must pass string of panel sides or dictionary mapping axes numbers to sides.')
    if not isinstance(innerpanels, (dict, str)):
        raise ValueError('Must pass string of panel sides or dictionary mapping axes numbers to sides.')
    innerpanels = axes_dict(innerpanels or '', False)
    innercolorbars = axes_dict(innercolorbars or '', False)
    if innercolorbars_kw:
        innerpanels_kw = innercolorbars_kw
    innerpanels_kw = axes_dict(innerpanels_kw, True)
    for kw in innerpanels_kw.values():
        kw['whichpanels'] = ''
    for num,which in innerpanels.items():
        innerpanels_kw[num]['whichpanels'] += translate(which)
    for num,which in innercolorbars.items():
        which = translate(which)
        if which:
            innerpanels_kw[num]['whichpanels'] += which
            if re.search('[bt]', which):
                kwargs['hspace'] = _default(kwargs.get('hspace',None), rc['gridspec.xlab'])
                innerpanels_kw[num]['sharex_panels'] = False
                innerpanels_kw[num]['hwidth'] = _default(innerpanels_kw[num].get('hwidth', None), rc['gridspec.cbar'])
                innerpanels_kw[num]['hspace'] = _default(innerpanels_kw[num].get('hspace', None), rc['gridspec.xlab'])
            if re.search('[lr]', which):
                kwargs['wspace'] = _default(kwargs.get('wspace',None), rc['gridspec.ylab'])
                innerpanels_kw[num]['sharey_panels'] = False
                innerpanels_kw[num]['wwidth'] = _default(innerpanels_kw[num].get('wwidth', None), rc['gridspec.cbar'])
                if 'l' in which and 'r' in which:
                    default = (rc['gridspec.ylab'], rc['gridspec.nolab'])
                elif 'l' in which:
                    default = rc['gridspec.ylab']
                else:
                    default = rc['gridspec.nolab']
                innerpanels_kw[num]['wspace'] = _default(innerpanels_kw[num].get('wspace', None), default)

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    figsize, offset, subplots_kw, gridspec_kw = \
            _subplots_kwargs(nrows, ncols, **kwargs)
    row_offset, col_offset = offset
    gs = gridspec.FlexibleGridSpec(**gridspec_kw)
    fig = plt.figure(figsize=figsize,
        gridspec=gs,
        subplots_kw=subplots_kw,
        auto_adjust=auto_adjust,
        pad=pad,
        rcreset=rcreset,
        silent=silent,
        FigureClass=FigureBase,
        )

    #--------------------------------------------------------------------------
    # Manage shared axes/axes with spanning labels
    #--------------------------------------------------------------------------
    # Get some axes properties
    # Note that these locations should be **sorted** by axes id
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty
    yrange = row_offset + np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = col_offset + np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    # asdfas
    # xmin   = np.array([xy[0].min() for xy in axes_ids]) # unused
    # ymax   = np.array([xy[1].max() for xy in axes_ids])

    # Shared axes: generate list of base axes-dependent axes pairs
    # That is, find where the minimum-maximum gridspec extent in 'x' for a
    # given axes matches the minimum-maximum gridspec extent for a base axes
    xgroups_base, xgroups_sorted, xgroups, grouped = [], [], [], []
    if sharex:
        for i in range(num_axes): # axes now have pseudo-numbers from 0 to num_axes-1
            matches       = (xrange[i,:]==xrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0] # gives ID number of matching_axes, from 0 to num_axes-1
            if i not in grouped and matching_axes.size>1:
                # Find all axes that have the same gridspec 'x' extents
                xgroups      += [matching_axes]
                # Get bottom-most axis with shared x; should be single number
                # xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
                xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
                # Sorted group
                xgroups_sorted += [matching_axes[np.argsort(yrange[matching_axes,1])[::-1]]] # bottom-most axes is first
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    ygroups_base, ygroups_sorted, ygroups, grouped = [], [], [], []
    if sharey:
        for i in range(num_axes):
            matches       = (yrange[i,:]==yrange).all(axis=1) # *broadcasting rules apply here*
            matching_axes = np.where(matches)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups      += [matching_axes]
                ygroups_base += [matching_axes[np.argmin(xrange[matching_axes,0])]] # left-most axis with shared y, for matching_axes
                ygroups_sorted += [matching_axes[np.argsort(xrange[matching_axes,0])]] # left-most axis is first
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already

    #--------------------------------------------------------------------------
    # Draw axes
    # TODO: Need to configure to automatically determine 'base' axes based on
    # what has already been drawn. Not critical but would be nice.
    # TODO: Need to do something similar for the spanning axes. Also will
    # allow label to be set on any of the axes, but when this happens, will
    # set the label on the 'base' spanning axes.
    #--------------------------------------------------------------------------
    # Base axes; to be shared with other axes as ._sharex, ._sharey attributes
    axs = num_axes*[None] # list of axes
    allgroups_base = []
    if sharex:
        allgroups_base += xgroups_base
    if sharey:
        allgroups_base += ygroups_base
    for i in allgroups_base:
        ax_kw = axes_kw[i]
        if axs[i] is not None: # already created
            continue
        if innerpanels_kw[i]['whichpanels']: # non-empty
            axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    spanx=spanx, spany=spany,
                    number=i+1, **ax_kw, **innerpanels_kw[i]) # main axes handle
        else:
            axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    spanx=spanx, spany=spany,
                    number=i+1, **ax_kw) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another. If so, get that
        # axes. Also do some error checking
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        ax_kw = axes_kw[i]
        if sharex:
            igroup = np.where([i in g for g in xgroups])[0]
            if igroup.size==1:
                sharex_ax = axs[xgroups_base[igroup[0]]]
                if sharex_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
        if sharey:
            igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
            if igroup.size==1:
                sharey_ax = axs[ygroups_base[igroup[0]]]
                if sharey_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')

        # Draw axes, and add to list
        if axs[i] is not None:
            # Axes is a *base* and has already been drawn, but might still
            # have shared axes (e.g. is bottom-axes of three-column plot
            # and we want it to share the leftmost y-axis)
            if sharex_ax is not None and axs[i] is not sharex_ax:
                axs[i]._sharex_setup(sharex_ax, sharex)
            if sharey_ax is not None and axs[i] is not sharey_ax:
                axs[i]._sharey_setup(sharey_ax, sharey)
        else:
            # Virgin axes; these are not an x base or a y base
            if innerpanels_kw[i]['whichpanels']: # non-empty
                axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1, spanx=spanx, spany=spany,
                        sharex_level=sharex, sharey_level=sharey,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw, **innerpanels_kw[i])
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        number=i+1, spanx=spanx, spany=spany,
                        sharex_level=sharex, sharey_level=sharey,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kw) # main axes can be a cartopy projection

    # Check that axes don't belong to multiple groups
    # This should be impossible unless my code is completely wrong...
    for ax in axs:
        for name,groups in zip(('sharex', 'sharey'), (xgroups, ygroups)):
            if sum(ax in group for group in xgroups)>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple {name} groups.')

    #--------------------------------------------------------------------------#
    # Create panel axes
    #--------------------------------------------------------------------------#
    def _paneladd(name, panels):
        if not panels:
            return
        axsp = []
        side = re.sub('^(.*)panel$', r'\1', name)
        for n in np.unique(panels).flat:
            offset = row_offset if side in ('left','right') else col_offset
            idx, = np.where(panels==n)
            idx = slice(offset + min(idx), offset + max(idx) + 1)
            if side=='right':
                subspec = gs[idx,-1]
            elif side=='left':
                subspec = gs[idx,0]
            elif side=='bottom':
                subspec = gs[-1,idx]
            axp = fig.add_subplot(subspec, panel_side=side, invisible=True, projection='panel')
            axsp += [axp]
        setattr(fig, name, axes_list(axsp))
    _paneladd('bottompanel', subplots_kw.bottompanels)
    _paneladd('rightpanel',  subplots_kw.rightpanels)
    _paneladd('leftpanel',   subplots_kw.leftpanels)

    # Return results
    # print('Figure setup complete.')
    return fig, axes_list(axs)

