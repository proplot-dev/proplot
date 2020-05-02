#!/usr/bin/env python3
"""
Default proplot configuration settings and validators for
`rc_configurator` assignments.
"""
# NOTE: Make sure to add to docs/configuration.rst when updating or adding
# new settings! Much of this script was adapted from seaborn; see:
# https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
from matplotlib import rcParams as _rc_params

# Dictionaries containing default settings
_rc_quick_defaults = {
    'abc': False,
    'align': False,
    'alpha': 1,
    'borders': False,
    'cmap': 'fire',
    'coast': False,
    'color': 'k',
    'cycle': 'colorblind',
    'facecolor': 'w',
    'fontname': 'sans-serif',
    'inlinefmt': 'retina',
    'geogrid': True,
    'grid': True,
    'gridminor': False,
    'gridratio': 0.5,
    'innerborders': False,
    'lakes': False,
    'land': False,
    'large': 10,
    'linewidth': 0.6,
    'lut': 256,
    'margin': 0.0,
    'ocean': False,
    'reso': 'lo',
    'rgbcycle': False,
    'rivers': False,
    'share': 3,
    'small': 9,
    'span': True,
    'tickdir': 'out',
    'ticklen': 4.0,
    'ticklenratio': 0.5,
    'tickminor': True,
    'tickpad': 2.0,
    'tickratio': 0.8,
    'tight': True,
}

_rc_added_defaults = {
    'abc.border': True,
    'abc.borderwidth': 1.5,
    'abc.color': 'k',
    'abc.loc': 'l',  # left side above the axes
    'abc.size': None,  # = large
    'abc.style': 'a',
    'abc.weight': 'bold',
    'axes.facealpha': None,  # if empty, depends on 'savefig.transparent'
    'axes.formatter.timerotation': 90,
    'axes.formatter.zerotrim': True,
    'axes.geogrid': True,
    'axes.gridminor': True,
    'borders.color': 'k',
    'borders.linewidth': 0.6,
    'bottomlabel.color': 'k',
    'bottomlabel.size': None,  # = large
    'bottomlabel.weight': 'bold',
    'coast.color': 'k',
    'coast.linewidth': 0.6,
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
    'geogrid.alpha': 0.5,
    'geogrid.color': 'k',
    'geogrid.labels': False,
    'geogrid.labelsize': None,  # = small
    'geogrid.latmax': 90,
    'geogrid.latstep': 20,
    'geogrid.linestyle': ':',
    'geogrid.linewidth': 1.0,
    'geogrid.lonstep': 30,
    'gridminor.alpha': None,  # = grid.alpha
    'gridminor.color': None,  # = grid.color
    'gridminor.linestyle': None,  # = grid.linewidth
    'gridminor.linewidth': None,  # = grid.linewidth x gridratio
    'image.edgefix': True,
    'image.levels': 11,
    'innerborders.color': 'k',
    'innerborders.linewidth': 0.6,
    'lakes.color': 'w',
    'land.color': 'k',
    'leftlabel.color': 'k',
    'leftlabel.size': None,  # = large
    'leftlabel.weight': 'bold',
    'ocean.color': 'w',
    'rightlabel.color': 'k',
    'rightlabel.size': None,  # = large
    'rightlabel.weight': 'bold',
    'rivers.color': 'k',
    'rivers.linewidth': 0.6,
    'subplots.axpad': '1em',
    'subplots.axwidth': '18em',
    'subplots.pad': '0.5em',
    'subplots.panelpad': '0.5em',
    'subplots.panelwidth': '4em',
    'suptitle.color': 'k',
    'suptitle.size': None,  # = large
    'suptitle.weight': 'bold',
    'tick.labelcolor': None,  # = color
    'tick.labelsize': None,  # = small
    'tick.labelweight': 'normal',
    'title.border': True,
    'title.borderwidth': 1.5,
    'title.color': 'k',
    'title.loc': 'c',  # centered above the axes
    'title.pad': 3.0,  # copy
    'title.size': None,  # = large
    'title.weight': 'normal',
    'toplabel.color': 'k',
    'toplabel.size': None,  # = large
    'toplabel.weight': 'bold',
}

_rc_param_defaults = {
    'axes.grid': True,
    'axes.labelpad': 3.0,
    'axes.titlepad': 3.0,
    'axes.titleweight': 'normal',
    'axes.xmargin': 0.0,
    'axes.ymargin': 0.0,
    'figure.autolayout': False,
    'figure.facecolor': '#f2f2f2',
    'figure.max_open_warning': 0,
    'figure.titleweight': 'bold',
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
    'grid.alpha': 0.1,
    'grid.color': 'k',
    'grid.linestyle': '-',
    'grid.linewidth': 0.6,
    'hatch.color': 'k',
    'hatch.linewidth': 0.6,
    'legend.borderaxespad': 0,
    'legend.borderpad': 0.5,
    'legend.columnspacing': 1.0,
    'legend.fancybox': False,
    'legend.framealpha': 0.8,
    'legend.frameon': True,
    'legend.handlelength': 1.5,
    'legend.handletextpad': 0.5,
    'legend.labelspacing': 0.5,
    'lines.linewidth': 1.3,
    'lines.markersize': 3.0,
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'savefig.bbox': 'standard',
    'savefig.directory': '',
    'savefig.dpi': 300,
    'savefig.facecolor': 'white',
    'savefig.format': 'pdf',
    'savefig.pad_inches': 0.0,
    'savefig.transparent': True,
    'text.usetex': False,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
}

# "Global" settings and the lower-level settings they change
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
    'large': (  # the 'large' fonts
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
    'ticklen': (
        'xtick.major.size', 'ytick.major.size'
    ),
    'tickdir': (
        'xtick.direction', 'ytick.direction'
    ),
    'labelpad': (
        'axes.labelpad',
    ),
    'titlepad': (
        'axes.titlepad',
    ),
    'tickpad': (
        'xtick.major.pad', 'xtick.minor.pad',
        'ytick.major.pad', 'ytick.minor.pad'
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
    for dict_ in (_rc_params, _rc_added_defaults, _rc_quick_defaults)
    for name in dict_.keys()
}

_rc_categories = {
    '.'.join(name.split('.')[:i + 1])
    for dict_ in (_rc_added_defaults, _rc_params)
    for name in dict_
    for i in range(len(name.split('.')) - 1)
}
