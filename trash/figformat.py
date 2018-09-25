def format(self, # axes=None, # optionally give list of axes
        # force=False, # override previous formatting
        mappable=None, handles=None, # for legends, colorbars
        titles=None, titlesinside=None, xlabels=None, ylabels=None, # common to want to apply all of these quickly
        filename=None, desktop=True, **kwargs): # the rest gets passed to the general format function
    """
    DELETED, BECAUSE THIS WORKFLOW IS KIND OF HARD TO FOLLOW AND UGLY; USER SHOULD JUST LOOP THROUGH EACH AXES
    OBJECT AND FORMAT THEM INDIVIDUALLY.
    Formatting multiple axes. Possible workflows...
        1) Format axes needing special treatment individually, then call this to format the rest.
        2) Format lists of axes explicitly with this function.
    ...doesn't really make sense to supply it with axes I think.
    """
    # Format list of axes
    # if axes is None:
    #     axes = self.axes
    #     # axes = self.axes_main
    axes = self.axes_main
    for i,a in enumerate(axes):
        # if not a.formatted or force:
        if titlesinside is not None:    kwargs.update(titleinside=titlesinside[i])
        if titles is not None:          kwargs.update(title=titles[i])
        if xlabels is not None:         kwargs.update(xlabel=xlabels[i])
        if ylabels is not None:         kwargs.update(ylabel=ylabels[i])
        a.format(**kwargs)

    # Format special axes
    if self.bottomlegend is not None and handles is not None:
        kwargs.update(handles=handles)
        self.bottomlegend.format(**kwargs)
    if self.rightcolorbar is not None and mappable is not None:
        kwargs.update(mappable=mappable)
        self.rightcolorbar.format(**kwargs)
    if self.bottomcolorbar is not None and mappable is not None:
        kwargs.update(mappable=mappable)
        self.bottomcolorbar.format(**kwargs)

    # And save, optionally (into figures folder)
    if filename is not None:
        self.save(filename, desktop=desktop)
