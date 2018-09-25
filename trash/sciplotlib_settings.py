def settings(ctick_kw=dict(), tick_kw=dict(), tickminor_kw=dict(), # tick stuff, useing set_tick_params
    ticklabels_kw=dict(), label_kw=dict(), title_kw=dict(), abc_kw=dict(), # font artists
    cgrid_kw=dict(), grid_kw=dict(), gridminor_kw=dict(), # line artists
    boundary_kw=dict(), lonlatlines_kw=dict(), continents_kw=dict(), coastlines_kw=dict(), # mapping line artists
    legend_kw=dict(), colorbar_kw=dict(), # special, used on legend/colmorbbar initiation
    xscale_kw=dict(), yscale_kw=dict(),
    color='k', linewidth=0.7, transparent=False, # axes options, spread to multiple artists
    ):
    """
    Default kwargs that are passed to format. Sort of like mpl_rc.
    Takes as kwargs several dictionaries containing properties to be **added**/**overridden**.
    """
    # Override defaults for legend/colorbar, labels, ticklabels, title, grid/gridminor
    # ...axes spines, ticks (set with set_tick_params)
    kw = Dict() # simple dictionary
    # ...scaling
    kw.xscale = xscale_kw # pass to set_xscale (see documentation; some special properties draw ticks easily)
    kw.yscale = yscale_kw # as above
    # ...and assign properties (length in points, here)
    kw.tick      = dict(dict(length=4, direction='in', width=linewidth, color=color), **tick_kw)
    kw.tickminor = dict(dict(length=2, direction='in', width=linewidth, color=color), **tickminor_kw) # dependent on above
    kw.ctick     = dict(dict(length=4, direction='out', width=linewidth, color=color), **ctick_kw)
    kw.spine      = dict(linewidth=linewidth, color=color) # spine keyword
    kw.box       = dict(linewidth=linewidth, edgecolor=color)
        # ...for outlines of colorbar, legend, cartopy frame, and cartopy features
    # ...grid lines (set on declaration)
    kw.grid      = dict(dict(linestyle='-', linewidth=0.8, color=color, alpha=0.1), **grid_kw)
    kw.gridminor = dict(dict(linestyle=':', linewidth=0.8, color=color, alpha=0.05), **gridminor_kw)
    # ...geographic (set on declaration)
    kw.coastlines  = dict(dict(linewidth=linewidth, color='#888888'), **coastlines_kw)
    kw.continents  = dict(dict(color='#CCCCCC'), **continents_kw)
    kw.lonlatlines = dict(dict(linestyle=':', linewidth=linewidth, color=color, alpha=0.2, dashes=(linewidth,linewidth)), **lonlatlines_kw)
    kw.boundary    = dict(dict(linewidth=2*linewidth, color=color), **boundary_kw)
        # boundary should be thick... for some reason looks skinnier than parallels unless multiply by 2
        # even though (checking after the fact) reveals the linewidths are held
    # ...text (set with plt.setp, set_label_text, or on creation)
    kw.ticklabels = dict(dict(size=8, weight='normal', color=color), **ticklabels_kw)
    kw.label      = dict(dict(size=8, weight='normal', color=color), **label_kw)
    kw.title      = dict(dict(size=10, weight='normal', color=color), **title_kw)
    kw.abc        = dict(dict(size=10, weight='bold', color=color), **abc_kw)
    # ...special (set on creation)
    kw.legend   = dict(dict(framealpha=0.6, fancybox=False, frameon=True, fontsize=8), **legend_kw)
    kw.colorbar = dict(dict(extend='both', spacing='uniform'), **colorbar_kw) # use cgrid, for ease
    kw.cgrid    = dict(dict(color=color, linewidth=linewidth), **cgrid_kw)
    return kw # returns the settings object

