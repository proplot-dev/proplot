#!/usr/bin/env python3
import matplotlib.gridspec as mgridspec

#------------------------------------------------------------------------------#
# Magical hacktastic class that delays drawing until user attempts looking up
# a method or something; useful when we want to fill an allocated space with
# something other than a typical axes (i.e. colorbar or legend)
#------------------------------------------------------------------------------#
# Helper methods because we are doing some seriously dubious fuckery here
get = object.__getattribute__
set_ = object.__setattr__

class Delayed(object):
    """
    Notes
    -----
    For inheriting from instantiated object see: https://stackoverflow.com/a/33468799/4970632
    but this only works for copying attributes, not methods.
    Goal of this is to have generalized panel class suitable for both 'inner'
    and 'outer' panels.
    """
    def __init__(self, figure, subspecs, instantiate=False, n=1, **kwargs):
        # Assignments
        npanel = n
        set_(self, '_axs',        [None]*n)
        set_(self, '_naxes',      n)
        set_(self, '_index',      0)
        set_(self, '_figure',     figure)
        set_(self, '_kwargs',     kwargs) # for instantiation
        set_(self, '_subspecs',   subspecs)
        set_(self, '_attributes', [{}]*n)
        # Optionally instantiate all the axes right away
        # Do this by setting the index each time
        if instantiate:
            get(self, 'instantiate')(**kwargs)

    def __setattr__(self, attr, value):
        index = get(self, '_index')
        get(self, '_attributes')[index].update({attr:value})

    def __getattribute__(self, attr, *args):
        """
        When method is first invoked, ensure axes are drawn. Depending on
        which property user accesses first, axes can be allocated as:
          * A legend box (just dead space for storing a legend)
          * A colorbar (optionally shrunken relative to subplotspec length)
          * A normal axes (for plotting and stuff)
        """
        # Instantiate axes
        init  = False
        index = get(self, '_index')
        axs   = get(self, '_axs')
        if not axs[index]: # use super-class method
            init = True
            get(self, 'instantiate')(attr)

        # Now that axes are initiated, get the relevant attribute
        # If lookup fails, delete the axes before raising the error
        if attr in ('colorbar','legend'):
            obj = self
        else:
            obj = axs[index]
        try:
            return get(obj, attr, *args)
        except AttributeError:
            if init:
                axs[index].set_visible(False)
                axs[index] = None # remove reference
            raise

    def instantiate(self, *args):
        """
        Function for instantiating axes belonging to this panel.
        Will read from the 'n' attribute to figure out how many panels to draw.
        """
        # This function is invoked from within __getattribute__
        # so need to make sure we don't trigger infinite loops
        axs    = get(self, '_axs')
        index  = get(self, '_index')
        if axs[index]:
            return axs[index]
        figure   = get(self, '_figure')
        subspecs = get(self, '_subspecs')

        # Instantiate
        kw = get(self, '_kwargs')
        projection = kw.pop('projection', 'xy')
        axs[index] = figure.add_subplot(subspecs[index], projection=projection, **kw)

        # Apply attributes to axes
        for attr,value in get(self, '_attributes')[index].items():
            ax = axs[index]
            setattr(ax, attr, value)
            if attr=='_sharex':
                for t in ax.xaxis.get_ticklabels():
                    t.set_visible(False)
                ax.xaxis.label.set_visible(False)
            if attr=='_sharey':
                for t in ax.yaxis.get_ticklabels():
                    t.set_visible(False)
                ax.yaxis.label.set_visible(False)

        return axs[index]

class DelayedAxes(Delayed):
    """
    As above, but for generalized axes.
    """
    def __init__(self, figure, subspec, **kwargs):
        # Assignments
        kwargs.update({'n':1}) # cannot have multiple axes here
        super().__init__(figure, [subspec], **kwargs)

    def legend(self, *args, **kwargs):
        # Dummy method
        index = get(self, '_index')
        axs   = get(self, '_axs')
        return axs[index].legend(*args, **kwargs)

class DelayedPanel(Delayed):
    """
    For panels, destined to be colorbars, axes, or legends. Behavior depends
    on which attribute you access first:
        * Try to access colorbar or legend, and this is a default
            axes.Axes instance.
        * Try to access anything else, and this is a normal XYAxes
            like everything else.
    Why do this? This way user can do whatever they want with panels, no
    extra kwargs necessary during subplots() call.
    """
    def __init__(self, figure, subspec, side, n=1, length=1, **kwargs):
        # Optionally shrink in the lengthwise direction and stack multiple
        # panels in the crosswise dimension
        if side in ('bottom','top'):
            gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows        = n,
                    ncols        = 3,
                    wspace       = 0,
                    hspace       = 0,
                    subplot_spec = subspec,
                    width_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            subspecs = [gridspec[i,1] for i in range(n)]
        elif side in ('left','right'):
            gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows         = 3,
                    ncols         = n,
                    hspace        = 0,
                    wspace        = 0,
                    subplot_spec  = subspec,
                    height_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            subspecs = [gridspec[1,i] for i in range(n)]
        else:
            raise ValueError(f'Invalid panel side "{side}".')

        # Call parent method
        set_(self, '_side', side)
        super().__init__(figure, subspecs, **kwargs)

    def __getitem__(self, key):
        naxes = get(self, '_naxes')
        if key<0:
            key = naxes + key # e.g. -1
        if key<0 or key>naxes-1:
            raise ValueError(f'Key {key} invalid, only {naxes} panels present.')
        self._index = key
        return self

    def instantiate(self, *args):
        # Override projection in some cases
        attr = None if len(args)==0 else args[0]
        if attr in ('colorbar','legend'):
            get(self, '_kwargs')['projection'] = None
        return super().instantiate(*args)

    def legend(self, handles, **kwargs):
        # Allocate invisible axes for drawing legend.
        # Returns the axes and the output of legend_factory().
        index = get(self, '_index')
        axs   = get(self, '_axs') # get axes
        ax    = axs[index]
        for s in ax.spines.values():
            s.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.patch.set_alpha(0)
        kwnew = {'frameon':False, 'loc':'center',
                 'borderaxespad':0,
                 'bbox_transform':ax.transAxes}
        kwnew.update(kwargs)
        return ax, legend_factory(ax, handles, **kwnew)

    def colorbar(self, *args, inner=False, outer=False, length=1, **kwargs):
        # Allocate axes for drawing colorbar.
        # Returns the axes and the output of colorbar_factory().
        index  = get(self, '_index')
        axs    = get(self, '_axs') # get axes
        ax     = axs[index]
        npanel = len(axs)
        side   = get(self, '_side')
        if side in ('bottom','top'):
            ticklocation = 'bottom' if index==npanel-1 else 'top'
            orientation  = 'horizontal'
        elif side in ('left','right'):
            ticklocation = 'right' if index==npanel-1 else 'right'
            orientation = 'vertical'
        kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})
        return ax, colorbar_factory(ax, *args, **kwargs)

