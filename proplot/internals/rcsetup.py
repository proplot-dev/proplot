#!/usr/bin/env python3
"""
Default proplot configuration settings and validators for
`rc_configurator` assignments.
"""
# NOTE: Make sure to add to docs/configuration.rst when updating or adding settings!
from matplotlib import rcParamsDefault as _rc_matplotlib_default


# Deprecated settings
# TODO: Add 'openrgb' setting for renaming shorthand color names to
_rc_removed = {'rgbcycle': '0.6'}  # {key: version} dictionary
_rc_renamed = {}  # {old_key: (new_key, version)} dictionary

# Quick settings
# TODO: This is currently a hodgepodge of "real" settings and meta-settings.
# Should give the real settings "real" longer names in _rc_added_default. For
# example 'land' should be 'land.draw' or 'land.on'.
# TODO: More consistent behavior for how format() handles rc params. Currently
# some things are only keyword arguments while others are actual settings, but not
# obvious which is which. For example xticklen and yticklen should be quick settings.
_rc_quick_default = {
    'abc': False,
    'align': False,
    'alpha': 1,
    'borders': False,
    'cmap': 'fire',
    'coast': False,
    'color': 'black',
    'cycle': 'colorblind',
    'facecolor': 'w',
    'fontname': 'sans-serif',
    'inlinefmt': 'retina',
    'geogrid': True,
    'grid': True,
    'gridminor': False,
    'innerborders': False,
    'lakes': False,
    'land': False,
    'large': 10,
    'linewidth': 0.8,
    'lut': 256,
    'margin': 0.05,
    'negcolor': 'blue7',
    'ocean': False,
    'poscolor': 'red7',
    'reso': 'med',
    'rivers': False,
    'share': 3,
    'small': 9,
    'span': True,
    'style': None,
    'tickminor': True,
    'tight': True,
}

# ProPlot pseudo-settings
# TODO: Implement these as bonafide matplotlib settings by subclassing
# matplotlib's RcParams and adding new validators. Quick settings should
# be implemented under __getitem__.
_rc_added_default = {
    'abc.border': True,
    'abc.borderwidth': 1.5,
    'abc.color': 'black',
    'abc.loc': 'l',  # left side above the axes
    'abc.size': None,  # = large
    'abc.style': 'a',
    'abc.weight': 'bold',
    'axes.facealpha': None,  # if empty, depends on 'savefig.transparent'
    'axes.formatter.timerotation': 90,
    'axes.formatter.zerotrim': True,
    'axes.geogrid': True,
    'axes.gridminor': True,
    'borders.color': 'black',
    'borders.linewidth': 0.8,
    'bottomlabel.color': 'black',
    'bottomlabel.size': None,  # = large
    'bottomlabel.weight': 'bold',
    'coast.color': 'black',
    'coast.linewidth': 0.8,
    'colorbar.extend': '1.3em',
    'colorbar.framealpha': 0.8,
    'colorbar.frameon': True,
    'colorbar.grid': False,
    'colorbar.insetextend': '1em',
    'colorbar.insetlength': '8em',
    'colorbar.insetpad': '0.5em',
    'colorbar.insetwidth': '1.2em',
    'colorbar.length': 1,
    'colorbar.loc': 'right',
    'colorbar.width': '1.5em',
    'geoaxes.edgecolor': None,  # = color
    'geoaxes.facealpha': None,  # = alpha
    'geoaxes.facecolor': None,  # = facecolor
    'geoaxes.linewidth': None,  # = linewidth
    'geogrid.alpha': 0.25,
    'geogrid.axisbelow': False,  # like axes.axisbelow
    'geogrid.color': 'black',
    'geogrid.labelpad': 5,  # use cartopy default
    'geogrid.labels': False,
    'geogrid.labelsize': None,  # = small
    'geogrid.latmax': 90,
    'geogrid.latstep': 20,
    'geogrid.linestyle': '-',
    'geogrid.linewidth': 0.8,
    'geogrid.lonstep': 30,
    'geogrid.rotatelabels': True,  # False limits projections where labels are available
    'grid.ratio': 0.5,
    'gridminor.alpha': None,  # = grid.alpha
    'gridminor.color': None,  # = grid.color
    'gridminor.linestyle': None,  # = grid.linewidth
    'gridminor.linewidth': None,  # = grid.linewidth x gridratio
    'image.edgefix': True,
    'image.levels': 11,
    'innerborders.color': 'black',
    'innerborders.linewidth': 0.8,
    'lakes.color': 'w',
    'land.color': 'black',
    'leftlabel.color': 'black',
    'leftlabel.size': None,  # = large
    'leftlabel.weight': 'bold',
    'ocean.color': 'w',
    'rightlabel.color': 'black',
    'rightlabel.size': None,  # = large
    'rightlabel.weight': 'bold',
    'rivers.color': 'black',
    'rivers.linewidth': 0.8,
    'subplots.axpad': '1em',
    'subplots.axwidth': '18em',
    'subplots.pad': '0.5em',
    'subplots.panelpad': '0.5em',
    'subplots.panelwidth': '4em',
    'suptitle.color': 'black',
    'suptitle.size': None,  # = large
    'suptitle.weight': 'bold',
    'tick.dir': 'out',
    'tick.labelcolor': None,  # = color
    'tick.labelpad': 2.0,
    'tick.labelsize': None,  # = small
    'tick.labelweight': 'normal',
    'tick.len': 4.0,
    'tick.lenratio': 0.5,
    'tick.pad': 2.0,
    'tick.ratio': 0.8,
    'title.border': True,
    'title.borderwidth': 1.5,
    'title.color': 'black',
    'title.loc': 'c',  # centered above the axes
    'title.pad': 3.0,  # copy
    'title.size': None,  # = large
    'title.weight': 'normal',
    'toplabel.color': 'black',
    'toplabel.size': None,  # = large
    'toplabel.weight': 'bold',
}

# ProPlot overrides of matplotlib default style
# NOTE: Hard to say what best value for 'margin' is. 0 is bad for bar plots and scatter
# plots, 0.05 is good for line plot in y direction but not x direction.
# NOTE: Settings bounds is same as setting limits, except bounds get
# overridden after next plot: https://stackoverflow.com/a/11467349/4970632
# NOTE: Some of these parameters are the same as matplotlib defaults but want
# to enforce some critical settings on top of user or system matplotlibrc files
_rc_params_default = {
    'axes.grid': True,  # enable lightweight transparent grid by default
    'axes.labelpad': 3.0,  # more compact
    'axes.titlepad': 3.0,  # more compact
    'axes.titleweight': 'normal',
    'figure.autolayout': False,
    'figure.facecolor': '#f2f2f2',  # similar to MATLAB interface
    'figure.titleweight': 'bold',  # differentiate from axes titles
    'font.serif': (
        'TeX Gyre Schola',  # Century lookalike
        'TeX Gyre Bonum',  # Bookman lookalike
        'TeX Gyre Termes',  # Times New Roman lookalike
        'TeX Gyre Pagella',  # Palatino lookalike
        'DejaVu Serif',
        'Bitstream Vera Serif',
        'Computer Modern Roman',
        'Bookman',
        'Century Schoolbook L',
        'Charter',
        'ITC Bookman',
        'New Century Schoolbook',
        'Nimbus Roman No9 L',
        'Palatino',
        'Times New Roman',
        'Times',
        'Utopia',
        'serif'
    ),
    'font.sans-serif': (
        'TeX Gyre Heros',  # Helvetica lookalike
        'DejaVu Sans',
        'Bitstream Vera Sans',
        'Computer Modern Sans Serif',
        'Arial',
        'Avenir',
        'Fira Math',
        'Frutiger',
        'Geneva',
        'Gill Sans',
        'Helvetica',
        'Lucid',
        'Lucida Grande',
        'Myriad Pro',
        'Noto Sans',
        'Roboto',
        'Source Sans Pro',
        'Tahoma',
        'Trebuchet MS',
        'Ubuntu',
        'Univers',
        'Verdana',
        'sans-serif'
    ),
    'font.monospace': (
        'TeX Gyre Cursor',  # Courier lookalike
        'DejaVu Sans Mono',
        'Bitstream Vera Sans Mono',
        'Computer Modern Typewriter',
        'Andale Mono',
        'Courier New',
        'Courier',
        'Fixed',
        'Nimbus Mono L',
        'Terminal',
        'monospace'
    ),
    'font.cursive': (
        'TeX Gyre Chorus',  # Chancery lookalike
        'Apple Chancery',
        'Felipa',
        'Sand',
        'Script MT',
        'Textile',
        'Zapf Chancery',
        'cursive'
    ),
    'font.fantasy': (
        'TeX Gyre Adventor',  # Avant Garde lookalike
        'Avant Garde',
        'Charcoal',
        'Chicago',
        'Comic Sans MS',
        'Futura',
        'Humor Sans',
        'Impact',
        'Optima',
        'Western',
        'xkcd',
        'fantasy'
    ),
    'grid.alpha': 0.1,  # lightweight unobtrusive gridlines
    'grid.color': 'black',  # lightweight unobtrusive gridlines
    'grid.linestyle': '-',
    'grid.linewidth': 0.8,
    'hatch.color': 'black',
    'hatch.linewidth': 0.8,
    'lines.linestyle': '-',
    'lines.linewidth': 1.5,
    'lines.markersize': 6.0,
    'legend.borderaxespad': 0,  # looks sleeker flush against edge
    'legend.borderpad': 0.5,  # a bit more space
    'legend.columnspacing': 1.5,  # more compact
    'legend.fancybox': False,  # looks modern without curvy box
    'legend.handletextpad': 0.5,
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'savefig.bbox': 'standard',  # use custom tight layout
    'savefig.directory': '',  # current directory
    'savefig.dpi': 300,  # low dpi to improve performance, high dpi when it matters
    'savefig.facecolor': 'white',  # different from figure.facecolor
    'savefig.format': 'pdf',  # most users use bitmap when vector graphics are better
    'savefig.transparent': True,
    'xtick.minor.visible': True,  # enable minor ticks by default
    'ytick.minor.visible': True,
}

# "Global" settings and the lower-level settings they change
# NOTE: Do not link title.color to axes.titlecolor because the latter
# can have value 'auto' which is not handled in format() right now,
# and because setting was only introduced in version 3.2.
_rc_children = {
    'cmap': (
        'image.cmap',
    ),
    'lut': (
        'image.lut',
    ),
    'alpha': (  # this is a custom setting
        'axes.facealpha', 'geoaxes.facealpha',
    ),
    'facecolor': (
        'axes.facecolor', 'geoaxes.facecolor'
    ),
    'fontname': (
        'font.family',
    ),
    'color': (  # change the 'color' of an axes
        'axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor',
        'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'
    ),
    'small': (  # the 'small' fonts
        'font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'
    ),
    'font.small': (
        'font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'
    ),
    'large': (  # the 'large' fonts
        'abc.size', 'figure.titlesize',
        'axes.titlesize', 'suptitle.size', 'title.size',
        'leftlabel.size', 'toplabel.size',
        'rightlabel.size', 'bottomlabel.size'
    ),
    'font.large': (
        'abc.size', 'figure.titlesize',
        'axes.titlesize', 'suptitle.size', 'title.size',
        'leftlabel.size', 'toplabel.size',
        'rightlabel.size', 'bottomlabel.size'
    ),
    'linewidth': (
        'axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth',
        'xtick.major.width', 'ytick.major.width'
    ),
    'margin': (
        'axes.xmargin', 'axes.ymargin'
    ),
    'grid': (
        'axes.grid',
    ),
    'gridminor': (
        'axes.gridminor',
    ),
    'geogrid': (
        'axes.geogrid',
    ),
    'tick.len': (
        'xtick.major.size', 'ytick.major.size'
    ),
    'tick.dir': (
        'xtick.direction', 'ytick.direction'
    ),
    'tick.labelpad': (
        'axes.labelpad',
    ),
    'tick.pad': (
        'xtick.major.pad', 'xtick.minor.pad',
        'ytick.major.pad', 'ytick.minor.pad'
    ),
    'title.pad': (
        'axes.titlepad',
    ),
    'grid.color': (
        'gridminor.color',
    ),
    'grid.linewidth': (
        'gridminor.linewidth',
    ),
    'grid.linestyle': (
        'gridminor.linestyle',
    ),
    'grid.alpha': (
        'gridminor.alpha',
    ),
}

# Settings without dots
_rc_nodots = {
    name.replace('.', ''): name
    for dict_ in (_rc_matplotlib_default, _rc_added_default, _rc_quick_default)
    for name in dict_.keys()
}

_rc_categories = {
    '.'.join(name.split('.')[:i + 1])
    for dict_ in (_rc_added_default, _rc_matplotlib_default)
    for name in dict_
    for i in range(len(name.split('.')) - 1)
}


# Helper functions
def _get_default_param(key):
    """
    Get the default parameter from one of three places. This is
    used for the :rc: role when compiling docs.
    """
    sentinel = object()
    for dict_ in (
        _rc_quick_default,
        _rc_added_default,
        _rc_params_default,
        _rc_matplotlib_default,
    ):
        value = dict_.get(key, sentinel)
        if value is not sentinel:
            return value
    raise KeyError(f'Invalid key {key!r}.')
