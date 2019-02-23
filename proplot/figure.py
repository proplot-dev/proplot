#!/usr/bin/env python3
"""
Figure subclass.
"""
import os
import numpy as np
import matplotlib.figure as mfigure
from .rcmod import rc
from .utils import _dot_dict, _default, _timer, _counter, units, ic
from . import gridspec, axes
# Aliases for panel names
_aliases = {
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'tpanel': 'toppanel',
    'lpanel': 'leftpanel'
    }

class Figure(mfigure.Figure):
    # Subclass adding some super cool features
    def __init__(self, figsize,
            gridspec=None, subplots_kw=None,
            rcreset=True,  auto_adjust=True, pad=0.1,
            **kwargs):
        """
        Matplotlib figure with some pizzazz.

        Parameters
        ----------
        figsize : (float or str, float or str)
            Figure size (width, height). If numeric, units are inches.
            Otherwise can specify alternative units.
        subplots_kw : None, dict-like
            Container of the keyword arguments used to initialize
            subplots.

        Other parameters
        ----------------
        rcreset : bool, optional
            When figure is drawn, reset rc settings to defaults?
        auto_adjust : bool, optional
            When figure is drawn, trim the gridspec edges without messing
            up axes aspect ratios and internal spacing?
        **kwargs : dict-like
            Passed to `matplotlib.figure.Figure` initializer.
        """
        # Initialize figure with some custom attributes.
        # Whether to reset rcParams wheenver a figure is drawn (e.g. after
        # ipython notebook finishes executing)
        self._rcreset  = rcreset
        self._smart_pad = units(pad)
        self._smart_tight = auto_adjust # note name _tight already taken!
        self._smart_tight_init = True # is figure in its initial state?
        self._span_labels = [] # add axis instances to this, and label position will be updated
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

    def __getattribute__(self, attr, *args):
        # Get attribute, but offer some convenient aliases
        attr = _aliases.get(attr, attr)
        return super().__getattribute__(attr, *args)

    def _rowlabels(self, labels, **kwargs):
        # Assign rowlabels
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
        # Assign collabels
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

    def _suptitle_setup(self, renderer=None, offset=False, **kwargs):
        # Intelligently determine supertitle position:
        # Determine x by the underlying gridspec structure, where main axes lie.
        left = self._subplots_kw.left
        right = self._subplots_kw.right
        if self.leftpanel:
            left += (self._subplots_kw.lwidth + self._subplots_kw.lspace)
        if self.rightpanel:
            right += (self._subplots_kw.rwidth + self._subplots_kw.rspace)
        xpos = left/self.width + 0.5*(self.width - left - right)/self.width

        if not offset or not kwargs.get('text', self._suptitle.get_text()):
            # Simple offset, not using the automatically determined
            # title position for guidance
            base = rc['axes.titlepad']/72 + self._gridspec.top*self.height
            ypos = base/self.height
            transform = self.transFigure
        else:
            # Figure out which title on the top-row axes will be offset the most
            # NOTE: Have to use private API to figure out whether axis has
            # tick labels or not! Seems to be no other way to do it.
            # See: https://matplotlib.org/_modules/matplotlib/axis.html#Axis.set_tick_params
            title_lev1, title_lev2, title_lev3 = None, None, None
            for ax in self.axes:
                # TODO: Need to ensure we do not test *bottom* axes panels
                if not isinstance(ax, axes.BaseAxes) or not ax._row_span[0]==0 or \
                    (isinstance(ax, axes.PanelAxes) and ax.panel_side=='bottom'):
                    continue
                title_lev1 = ax.title # always will be non-None
                if ((ax.title.get_text() and not ax._title_inside) or ax.collabel.get_text()):
                    title_lev2 = ax.title
                if ax.xaxis.get_ticks_position() == 'top':
                    test = 'label1On' not in ax.xaxis._major_tick_kw \
                        or ax.xaxis._major_tick_kw['label1On'] \
                        or ax.xaxis._major_tick_kw['label2On']
                    if test:
                        title_lev3 = ax.title

            # Hacky bugfixes:
            # 1) If no title, fill with spaces. Does nothing in most cases, but
            # if tick labels are on top, without this step matplotlib tight subplots
            # will not see the suptitle; now suptitle will just occupy empty title space.
            # 2) If title and tick labels on top, offset the suptitle and get
            # matplotlib to adjust tight_subplot by prepending newlines to title.
            # 3) Otherwise, offset suptitle, and matplotlib will recognize the
            # suptitle during tight_subplot adjustment.
            # ic(title_lev2, title_lev1, title_lev3)
            if not title_lev2: # no title present
                line = 0
                title = title_lev1
                if not title.get_text():
                # if not title.axes._title_inside:
                    title.set_text('\n ') # dummy spaces, so subplots adjust will work properly
            elif title_lev3: # upper axes present
                line = 1.0 # looks best empirically
                title = title_lev3
                text = title.get_text()
                title.set_text('\n\n' + text)
            else:
                line = 1.2 # default line spacing; see: https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
                title = title_lev1 # most common one

            # First idea: Create blended transform, end with newline
            # ypos = title.get_position()[1]
            # transform = mtransforms.blended_transform_factory(
            #         self.transFigure, title.get_transform())
            # text = kwargs.pop('text', self._suptitle.get_text())
            # if text[-1:] != '\n':
            #     text += '\n'
            # kwargs['text'] = text
            # New idea: Get the transformed position
            # NOTE: Seems draw() is called more than once, and the last times
            # are when title positions are appropriately offset.
            # NOTE: Default linespacing is 1.2; it has no get, only a setter; see
            # https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
            transform = title.get_transform() + self.transFigure.inverted()
            ypos = transform.transform(title.get_position())[1]
            line = line*(rc['axes.titlesize']/72)/self.height
            ypos = ypos + line
            transform = self.transFigure

        # Update settings
        self._suptitle.update({'position':(xpos, ypos),
            'transform':transform,
            'ha':'center', 'va':'bottom', **kwargs})

    def draw(self, renderer, *args, **kwargs):
        """
        Custom draw.
        """
        # Special: Figure out if other titles are present, and if not
        # bring suptitle close to center
        # ref_ax = None
        self._suptitle_setup(renderer, offset=True) # just applies the spacing
        self._auto_smart_tight_layout(renderer)
        # If rc settings have been changed, reset them when the figure is
        # displayed (usually means we have finished executing a notebook cell).
        if not rc._init and self._rcreset:
            print('Resetting rcparams.')
            rc.reset()
        return super().draw(renderer, *args, **kwargs)

    def panel_factory(self, subspec, whichpanels=None,
            hspace=None, wspace=None,
            hwidth=None, wwidth=None,
            sharex=None, sharey=None, # external sharing
            sharex_level=3, sharey_level=3,
            sharex_panels=True, sharey_panels=True, # by default share main x/y axes with panel x/y axes
            **kwargs):
        """
        Create edge panels.
        """
        # Helper function for creating paneled axes.
        width, height = self.width, self.height
        translate = {'bottom':'b', 'top':'t', 'right':'r', 'left':'l'}
        whichpanels = translate.get(whichpanels, whichpanels)
        whichpanels = whichpanels or 'r'
        hspace = _default(hspace, 0.13) # teeny tiny space
        wspace = _default(wspace, 0.13)
        hwidth = _default(hwidth, 0.45) # default is panels for plotting stuff, not colorbars
        wwidth = _default(wwidth, 0.45)
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
        hwidth_ratios = [height-hwidth*(nrows-1)]*nrows
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

    def smart_tight_layout(self, renderer=None, pad=None):
        """
        Get arguments necessary passed to subplots() to create a tight figure
        bounding box without screwing aspect ratios, widths/heights, and such.
        """
        # Get bounding box that encompasses *all artists*, compare to bounding
        # box used for saving *figure*
        if pad is None:
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
        figsize, *_, gridspec_kw = gridspec._gridspec_kwargs(**subplots_kw)
        self._smart_tight_init = False
        self._gridspec.update(**gridspec_kw)
        self.set_size_inches(figsize)

        # Fix any spanning labels that we've added to _span_labels
        # These need figure-relative height coordinates
        for axis in self._span_labels:
            axis.axes._share_span_label(axis)

    def _auto_smart_tight_layout(self, renderer=None):
        # If we haven't already, compress edges
        if not self._smart_tight_init or not self._smart_tight:
            return
        # Cartopy sucks at labels! Bounding box identified will be wrong.
        # 1) If you used set_bounds to zoom into part of a cartopy projection,
        # this can erroneously identify invisible edges of map as being part of boundary
        # 2) If you have gridliner text labels, matplotlib won't detect them.
        if any(isinstance(ax, axes.CartopyAxes) for ax in self.axes):
            return
        # Adjust if none of not done already
        print('Adjusting gridspec.')
        self.smart_tight_layout(renderer)

    def save(self, filename, silent=False, auto_adjust=True, pad=0.1, **kwargs):
        """
        Custom save.
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
            kwargs['transparent'] = True
        # Finally, save
        self._suptitle_setup(offset=True) # just applies the spacing
        self._auto_smart_tight_layout()
        if not silent:
            print(f'Saving to "{filename}".')
        return super().savefig(os.path.expanduser(filename), **kwargs) # specify DPI for embedded raster objects

    def savefig(self, *args, **kwargs):
        """
        Custom save.
        """
        # Alias for save.
        return self.save(*args, **kwargs)

