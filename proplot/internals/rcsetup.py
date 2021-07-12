#!/usr/bin/env python3
"""
Default proplot configuration settings and validators for
`RcConfigurator` assignments.
"""
import numbers

import cycler
import numpy as np
from matplotlib import rcParamsDefault as _rc_matplotlib_default_full

from . import warnings

# Initial synced properties
# NOTE: Important that LINEWIDTH is less than matplotlib default of 0.8.
# In general want axes lines to look about as thick as text.
COLOR = 'black'
CMAP = 'fire'
CYCLE = 'colorblind'
CYCLIC = 'twilight'
DIVERGING = 'div'
FRAMEALPHA = 0.8  # legend and colorbar
FONTSIZE = 9.0
GRIDALPHA = 0.11
GRIDBELOW = 'line'
GRIDCOLOR = 'black'
GRIDRATIO = 0.5  # differentiated from major by half size reduction
GRIDSTYLE = '-'
LABELPAD = 4.0  # default is 4.0, previously was 3.0
LABELSIZE = 'medium'
LINEWIDTH = 0.6
MARGIN = 0.05
MATHTEXT = False
TICKDIR = 'out'
TICKLEN = 4.0
TICKLENRATIO = 0.5  # differentiated from major by half length reduction
TICKMINOR = True
TICKPAD = 2.0
TICKRATIO = 0.8  # very slight width reduction
TITLEPAD = 5.0  # default is 6.0, previously was 3.0
TITLESIZE = 'med-large'
ZLINES = 2  # default zorder for lines
ZPATCHES = 1

# Deprecated settings
# TODO: Add 'openrgb' setting for renaming shorthand color names to open-color colors
_rc_removed = {  # {key: (alternative, version)} dictionary
    'rgbcycle': ('', '0.6'),  # no alternative, we no longer offer this feature
    'geogrid.lonstep': ('Use ax.format(lonlocator=N) instead.', '0.6'),
    'geogrid.latstep': ('Use ax.format(latlocator=N) instead.', '0.6'),
    'geogrid.latmax': ('Use ax.format(latmax=N) instead.', '0.6'),
}
_rc_renamed = {  # {old_key: (new_key, version)} dictionary
    'abc.format': ('abc.style', '0.5'),
    'align': ('subplots.align', '0.6'),
    'axes.facealpha': ('axes.alpha', '0.6'),
    'geoaxes.edgecolor': ('axes.edgecolor', '0.6'),
    'geoaxes.facealpha': ('axes.alpha', '0.6'),
    'geoaxes.facecolor': ('axes.facecolor', '0.6'),
    'geoaxes.linewidth': ('axes.linewidth', '0.6'),
    'geogrid.alpha': ('grid.alpha', '0.6'),
    'geogrid.color': ('grid.color', '0.6'),
    'geogrid.labels': ('grid.labels', '0.6'),
    'geogrid.labelpad': ('grid.pad', '0.6'),
    'geogrid.labelsize': ('grid.labelsize', '0.6'),
    'geogrid.linestyle': ('grid.linestyle', '0.6'),
    'geogrid.linewidth': ('grid.linewidth', '0.6'),
    'share': ('subplots.share', '0.6'),
    'small': ('text.labelsize', '0.6'),
    'large': ('text.titlesize', '0.6'),
    'span': ('subplots.span', '0.6'),
    'tight': ('subplots.tight', '0.6'),
    'tick.labelpad': ('tick.pad', '0.6'),
    'axes.formatter.timerotation': ('formatter.timerotation', '0.6'),
    'axes.formatter.zerotrim': ('formatter.zerotrim', '0.6'),
    'abovetop': ('title.above', '0.7'),
    'subplots.pad': ('subplots.outerpad', '0.7'),
    'subplots.axpad': ('subplots.innerpad', '0.7'),
    'subplots.axwidth': ('subplots.refwidth', '0.7'),
}

# ProPlot overrides of matplotlib default style
# NOTE: Hard to say what best value for 'margin' is. 0 is bad for bar plots and scatter
# plots, 0.05 is good for line plot in y direction but not x direction.
# WARNING: Critical to include every parameter here that can be changed by a
# "meta" setting so that _get_default_param returns the value imposed by *proplot*
# and so that "changed" settings detectd by RcConfigurator.save are correct.
_rc_matplotlib_default = {
    'axes.axisbelow': GRIDBELOW,
    'axes.formatter.use_mathtext': MATHTEXT,
    'axes.grid': True,  # enable lightweight transparent grid by default
    'axes.grid.which': 'major',
    'axes.labelpad': LABELPAD,  # more compact
    'axes.labelsize': LABELSIZE,
    'axes.labelweight': 'normal',
    'axes.linewidth': LINEWIDTH,
    'axes.titlepad': TITLEPAD,  # more compact
    'axes.titlesize': TITLESIZE,
    'axes.titleweight': 'normal',
    'axes.xmargin': MARGIN,
    'axes.ymargin': MARGIN,
    'figure.autolayout': False,
    'figure.dpi': 100,
    'figure.facecolor': '#f2f2f2',  # similar to MATLAB interface
    'figure.titlesize': TITLESIZE,
    'figure.titleweight': 'bold',  # differentiate from axes titles
    'font.serif': [  # NOTE: font lists passed to rcParams are lists, not tuples
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
    ],
    'font.sans-serif': [
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
    ],
    'font.monospace': [
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
    ],
    'font.cursive': [
        'TeX Gyre Chorus',  # Chancery lookalike
        'Apple Chancery',
        'Felipa',
        'Sand',
        'Script MT',
        'Textile',
        'Zapf Chancery',
        'cursive'
    ],
    'font.fantasy': [
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
    ],
    'font.size': FONTSIZE,
    'grid.alpha': GRIDALPHA,  # lightweight unobtrusive gridlines
    'grid.color': GRIDCOLOR,  # lightweight unobtrusive gridlines
    'grid.linestyle': GRIDSTYLE,
    'grid.linewidth': LINEWIDTH,
    'hatch.color': COLOR,
    'hatch.linewidth': LINEWIDTH,
    'image.cmap': CMAP,
    'lines.linestyle': '-',
    'lines.linewidth': 1.5,
    'lines.markersize': 6.0,
    'legend.borderaxespad': 0,  # looks sleeker flush against edge
    'legend.borderpad': 0.5,  # a bit more space
    'legend.columnspacing': 1.5,  # more compact
    'legend.fancybox': False,  # looks modern without curvy box
    'legend.fontsize': LABELSIZE,
    'legend.framealpha': FRAMEALPHA,
    'legend.handletextpad': 0.5,
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'patch.linewidth': LINEWIDTH,
    'savefig.bbox': None,  # use custom tight layout
    'savefig.directory': '',  # current directory
    'savefig.dpi': 1000,  # academic journal recommendations for raster line art
    'savefig.facecolor': 'white',  # different from figure.facecolor
    'savefig.format': 'pdf',  # most users use bitmap, but vector graphics are better
    'savefig.transparent': False,
    'xtick.direction': TICKDIR,
    'xtick.labelsize': LABELSIZE,
    'xtick.major.pad': TICKPAD,
    'xtick.major.size': TICKLEN,
    'xtick.major.width': LINEWIDTH,
    'xtick.minor.pad': TICKPAD,
    'xtick.minor.size': TICKLEN * TICKLENRATIO,
    'xtick.minor.visible': TICKMINOR,
    'xtick.minor.width': LINEWIDTH * TICKRATIO,
    'ytick.direction': TICKDIR,
    'ytick.labelsize': LABELSIZE,
    'ytick.major.pad': TICKPAD,
    'ytick.major.size': TICKLEN,
    'ytick.major.width': LINEWIDTH,
    'ytick.minor.pad': TICKPAD,
    'ytick.minor.size': TICKLEN * TICKLENRATIO,
    'ytick.minor.width': LINEWIDTH * TICKRATIO,
    'ytick.minor.visible': TICKMINOR,
}

# Proplot pseudo-settings
# TODO: More consistent behavior for how format() handles rc params. Currently
# some things are only keyword arguments while others are actual settings, but not
# obvious which is which. For example xticklen and yticklen should be quick settings.
# TODO: Implement these as bonafide matplotlib settings by subclassing matplotlib's
# RcParams and adding validators. Quick settings should be implemented in __getitem__.
# NOTE: Cannot have different a-b-c and title paddings because they are both controlled
# by matplotlib's _title_offset_trans transform and want to keep them aligned anyway.
_addendum_units = ' Interpreted by `~proplot.utils.units`.'
_addendum_fonts = (
    ' (see `this list of valid font sizes '
    '<https://matplotlib.org/stable/tutorials/text/text_props.html#default-font>`__).'
)
_rc_proplot = {
    # Stylesheet
    'style': (
        None,
        'The default matplotlib `stylesheet '
        '<https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html>`__ '  # noqa: E501
        'name. If ``None``, a custom proplot style is used. '
        "If ``'default'``, the default matplotlib style is used."
    ),

    # A-b-c labels
    'abc': (
        False,
        'Boolean, whether to draw a-b-c labels by default.'
    ),
    'abc.border': (
        True,
        'Boolean, indicates whether to draw a white border around a-b-c labels '
        'when :rcraw:`abc.loc` is inside the axes.'
    ),
    'abc.borderwidth': (
        1.5,
        'Width of the white border around a-b-c labels.'
    ),
    'abc.bbox': (
        False,
        'Boolean, whether to draw semi-transparent bounding boxes around a-b-c labels '
        'when :rcraw:`abc.loc` is inside the axes.'
    ),
    'abc.bboxcolor': (
        'w',
        'a-b-c label bounding box color.'
    ),
    'abc.bboxstyle': (
        'square',
        'a-b-c label bounding box style.'
    ),
    'abc.bboxalpha': (
        0.5,
        'a-b-c label bounding box opacity.'
    ),
    'abc.bboxpad': (
        None,
        'Padding for the a-b-c label bounding box. By default this is scaled '
        'to make the box flush against the subplot edge. ' + _addendum_units
    ),
    'abc.color': (
        'black',
        'a-b-c label color.'
    ),
    'abc.loc': (
        'l',  # left side above the axes
        'a-b-c label position. For options, see the :ref:`title location '
        'table <title_table>`.'
    ),
    'abc.size': (
        TITLESIZE,
        'a-b-c label font size.'
    ),
    'abc.style': (
        'a',
        'a-b-c label style. Must be string containing the character ``a`` '
        "or ``A``, for example ``'a.'`` or ``'(A)'``."
    ),
    'abc.titlepad': (
        '0.5em',
        'Padding used to displace the title and the a-b-c label when they are '
        'in the same location.'
    ),
    'abc.weight': (
        'bold',
        'a-b-c label font weight.'
    ),
    'autoformat': (
        True,
        'Whether to automatically apply labels from `pandas.Series`, '
        '`pandas.DataFrame`, and `xarray.DataArray` objects passed to '
        'plotting functions.'
    ),

    # Axes additions
    'alpha': (
        1.0,
        'The opacity of the background axes patch.'
    ),
    'axes.alpha': (
        1.0,
        'The opacity of the background axes patch.'
    ),
    'axes.titleabove': (
        True,
        'Boolean, indicates whether to move the title and a-b-c labels above any "top" '
        'panels above axes.'
    ),
    'formatter.timerotation': (
        90,
        'Float, indicates the default *x* axis tick label rotation '
        'for datetime tick labels.'
    ),
    'formatter.zerotrim': (
        True,
        'Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.'
    ),
    'formatter.limits': (
        (-5, 6),
        'Alias for :rcraw:`axes.formatter.limits`.'
    ),
    'formatter.use_locale': (
        False,
        'Alias for :rcraw:`axes.formatter.use_locale`.'
    ),
    'formatter.use_mathtext': (
        MATHTEXT,
        'Alias for :rcraw:`axes.formatter.use_mathtext`.'
    ),
    'formatter.min_exponent': (
        0,
        'Alias for :rcraw:`axes.formatter.min_exponent`.'
    ),
    'formatter.use_offset': (
        True,
        'Alias for :rcraw:`axes.formatter.useOffset`.'
    ),
    'formatter.offset_threshold': (
        4,
        'Alias for :rcraw:`axes.formatter.offset_threshold`.'
    ),

    # Special basemap settings
    'basemap': (
        False,
        'Boolean, toggles whether basemap is the default backend.'
    ),

    # Country borders
    'borders': (
        False,
        'Boolean, toggles country border lines on and off.'
    ),
    'borders.color': (
        'black',
        'Line color for country borders.'
    ),
    'borders.linewidth': (
        LINEWIDTH,
        'Line width for country borders.'
    ),
    'borders.zorder': (
        ZLINES,
        'Z-order for country border lines.'
    ),

    # Bottom subplot labels
    'bottomlabel.color': (
        'black',
        'Font color for column labels on the bottom of the figure.'
    ),
    'bottomlabel.size': (
        TITLESIZE,
        'Font size for column labels on the bottom of the figure.'
    ),
    'bottomlabel.weight': (
        'bold',
        'Font weight for column labels on the bottom of the figure.'
    ),
    'bottomlabel.pad': (
        '0.3em',
        'Padding between axes content and column labels on the bottom of the figure. '
        + _addendum_units
    ),

    # Special cartopy settings
    'cartopy.autoextent': (
        False,
        'If ``False`` (the default), cartopy projection extents are global by '
        'default and no longer automatically adjusted based on plotted content.'
    ),
    'cartopy.circular': (
        True,
        "If ``True`` (the default), polar cartopy projections like ``'npstere'`` and "
        "``'spstere'`` are bounded with circles rather than squares."
    ),

    # Coastlines
    'coast': (
        False,
        'Boolean, toggles coastline lines on and off.'
    ),
    'coast.color': (
        'black',
        'Line color for coast lines.'
    ),
    'coast.linewidth': (
        LINEWIDTH,
        'Line width for coast lines.'
    ),

    # Colorbars
    'colorbar.extend': (
        '1.3em',
        'Length of rectangular or triangular "extensions" for panel colorbars. '
        + _addendum_units
    ),
    'colorbar.framealpha': (
        FRAMEALPHA,
        'Opacity for inset colorbar frames.'
    ),
    'colorbar.frameon': (
        True,
        'Boolean, indicates whether to draw a frame behind inset colorbars.'
    ),
    'colorbar.grid': (
        False,
        'Boolean, indicates whether to draw borders between each level of the colorbar.'
    ),
    'colorbar.insetextend': (
        '1em',
        'Length of rectangular or triangular "extensions" for inset colorbars. '
        + _addendum_units
    ),
    'colorbar.insetlength': (
        '8em',
        'Length of inset colorbars. ' + _addendum_units
    ),
    'colorbar.insetpad': (
        '0.7em',
        'Padding between axes edge and inset colorbars. ' + _addendum_units
    ),
    'colorbar.insetwidth': (
        '1.2em',
        'Width of inset colorbars. ' + _addendum_units
    ),
    'colorbar.length': (
        1,
        'Length of outer colorbars.'
    ),
    'colorbar.loc': (
        'right',
        'Inset colorbar location. For options, see the :ref:`location table '
        '<colorbar_table>`.'
    ),
    'colorbar.width': (
        '1.5em',
        'Width of outer colorbars. ' + _addendum_units
    ),

    # Style shorthands
    'cmap': (
        CMAP,
        'The default sequential colormap.'
    ),
    'color': (
        COLOR,
        'The color of axis spines, tick marks, tick labels, and labels.'
    ),
    'cycle': (
        CYCLE,
        'The name of the color cycle used for plot elements like lines.',
    ),
    'facecolor': (
        'white',
        'The color of the background axes patch.'
    ),

    # Font settings
    'font.name': (
        'sans-serif',
        "Alias for :rcraw:`font.family`. The default is ``'sans-serif'``."
    ),

    # Gridlines
    'grid': (
        True,
        'Boolean, toggles major grid lines on and off.'
    ),
    'grid.below': (
        GRIDBELOW,  # like axes.axisbelow
        'Alias for :rcraw:`axes.axisbelow`. If ``False``, draw gridlines on top of '
        "everything. If ``True``, underneath everything. If ``'line'``, "
        'underneath patches only.'
    ),
    'grid.dmslabels': (
        True,
        'Boolean, indicates whether to use degrees-minutes-seconds rather than '
        'decimals for gridline labels on `~proplot.axes.CartopyAxes`.'
    ),
    'grid.pad': (
        5,
        'Padding between map boundary edge and longitude and '
        'latitude labels for `~proplot.axes.GeoAxes`. ' + _addendum_units
    ),
    'grid.labels': (
        False,
        'Boolean, indicates whether to label the longitude and latitude gridlines '
        'in `~proplot.axes.GeoAxes`.'
    ),
    'grid.labelsize': (
        LABELSIZE,
        'Font size for longitude and latitude gridline labels in '
        '`~proplot.axes.GeoAxes`.'
    ),
    'grid.labelweight': (
        'normal',
        'Font weight for longitude and latitude gridline labels in '
        '`~proplot.axes.GeoAxes`.'
    ),
    'grid.labelcolor': (
        COLOR,
        'Font color for longitude and latitude gridline labels in '
        '`~proplot.axes.GeoAxes`.'
    ),
    'grid.latinline': (
        False,
        'Whether to use inline labels for `~proplot.axes.CartopyAxes` '
        'latitude gridlines.'
    ),
    'grid.loninline': (
        False,
        'Whether to use inline labels for `~proplot.axes.CartopyAxes` '
        'longitude gridlines.'
    ),
    'grid.nsteps': (
        250,
        'Number of interpolation steps used to draw cartopy gridlines.'
    ),
    'grid.ratio': (
        GRIDRATIO,
        'Ratio of minor gridline width to major gridline width.'
    ),
    'grid.rotatelabels': (
        False,  # False limits projections where labels are available
        'Boolean, indicates whether to rotate longitude and latitude '
        '`~proplot.axes.CartopyAxes` gridline labels.'
    ),

    # Minor gridlines
    'gridminor': (
        False,
        'Boolean, toggles minor grid lines on and off.'
    ),
    'gridminor.alpha': (
        GRIDALPHA,
        'Minor gridline transparency.'
    ),
    'gridminor.color': (
        GRIDCOLOR,
        'Minor gridline color.'
    ),
    'gridminor.latstep': (
        10,
        'Latitude gridline interval for `~proplot.axes.GeoAxes` with global '
        'extent.'
    ),
    'gridminor.linestyle': (
        GRIDSTYLE,
        'Minor gridline style.'
    ),
    'gridminor.linewidth': (
        GRIDRATIO * LINEWIDTH,
        'Minor gridline width.'
    ),
    'gridminor.lonstep': (
        20,
        'Interval for `~proplot.axes.GeoAxes` longitude gridlines in degrees.'
    ),

    # Image additions
    'image.discrete': (
        None,
        'If ``True``, `~proplot.colors.DiscreteNorm` is used for every colormap plot. '
        'If ``False``, it is never used. If ``None``, it is used for all plot types '
        'except `imshow`, `matshow`, `spy`, `hexbin`, and `hist2d`.'
    ),
    'image.edgefix': (
        True,
        'Whether to fix the `white-lines-between-filled-contours '
        '<https://stackoverflow.com/q/8263769/4970632>`__ and '
        '`white-lines-between-pcolor-rectangles '
        '<https://stackoverflow.com/q/27092991/4970632>`__ issues.'
    ),
    'image.inbounds': (
        True,
        'If ``True`` and the *x* and *y* axis limits have been explicitly set, '
        'only in-bounds data is considered when determining default colormap limits.'
    ),
    'image.levels': (
        11,
        'Default number of `~proplot.colors.DiscreteNorm` levels for plotting '
        'commands that use colormaps.'
    ),

    # Backend stuff
    'inlinefmt': (
        'retina',
        'The inline backend figure format or list thereof. Valid formats include '
        "``'svg'``, ``'pdf'``, ``'retina'``, ``'png'``, and ``jpeg``."
    ),

    # Inner borders
    'innerborders': (
        False,
        'Boolean, toggles internal political border lines (e.g. states and provinces) '
        'on and off.'
    ),
    'innerborders.color': (
        'black',
        'Line color for internal political borders.'
    ),
    'innerborders.linewidth': (
        LINEWIDTH,
        'Line width for internal political borders.'
    ),
    'innerborders.zorder': (
        ZLINES,
        'Z-order for internal border lines.'
    ),

    # Lake patches
    'lakes': (
        False,
        'Boolean, toggles lake patches on and off.'
    ),
    'lakes.color': (
        'w',
        'Face color for lake patches.'
    ),
    'lakes.zorder': (
        ZPATCHES,
        'Z-order for lake patches.'
    ),

    # Land patches
    'land': (
        False,
        'Boolean, toggles land patches on and off.'
    ),
    'land.color': (
        'black',
        'Face color for land patches.'
    ),
    'land.zorder': (
        ZPATCHES,
        'Z-order for land patches.'
    ),

    # Left subplot labels
    'leftlabel.color': (
        'black',
        'Font color for row labels on the left-hand side.'
    ),
    'leftlabel.size': (
        TITLESIZE,
        'Font size for row labels on the left-hand side.'
    ),
    'leftlabel.weight': (
        'bold',
        'Font weight for row labels on the left-hand side.'
    ),
    'leftlabel.pad': (
        '0.6em',
        'Padding between axes content and row labels on the left-hand side. '
        + _addendum_units
    ),

    # Edge width bulk setting
    'linewidth': (
        LINEWIDTH,
        'Thickness of axes spines and major tick lines.'
    ),

    # Misc bulk settings
    'lut': (
        256,
        'The number of colors to put in the colormap lookup table.'
    ),
    'margin': (
        MARGIN,
        'The margin of space between axes edges and objects plotted inside the axes '
        'if ``xlim`` and ``ylim`` are unset.'
    ),

    # For negative positive patches
    'negcolor': (
        'blue7',
        'The color for negative bars and shaded areas when using ``negpos=True``. '
        'See also :rcraw:`poscolor`.'
    ),
    'poscolor': (
        'red7',
        'The color for positive bars and shaded areas when using ``negpos=True``. '
        'See also :rcraw:`negcolor`.'
    ),

    # Ocean patches
    'ocean': (
        False,
        'Boolean, toggles ocean patches on and off.'
    ),
    'ocean.color': (
        'w',
        'Face color for ocean patches.'
    ),
    'ocean.zorder': (
        ZPATCHES,
        'Z-order for ocean patches.'
    ),

    # Geographic resolution
    'reso': (
        'lo',
        'Resolution for `~proplot.axes.GeoAxes` geographic features. '
        "Must be one of ``'lo'``, ``'med'``, ``'hi'``, ``'x-hi'``, or ``'xx-hi'``."
    ),

    # Right subplot labels
    'rightlabel.color': (
        'black',
        'Font color for row labels on the right-hand side.'
    ),
    'rightlabel.size': (
        TITLESIZE,
        'Font size for row labels on the right-hand side.'
    ),
    'rightlabel.weight': (
        'bold',
        'Font weight for row labels on the right-hand side.'
    ),
    'rightlabel.pad': (
        '0.6em',
        'Padding between axes content and row labels on the right-hand side. '
        + _addendum_units
    ),

    # River lines
    'rivers': (
        False,
        'Boolean, toggles river lines on and off.'
    ),
    'rivers.color': (
        'black',
        'Line color for river lines.'
    ),
    'rivers.linewidth': (
        LINEWIDTH,
        'Line width for river lines.'
    ),
    'rivers.zorder': (
        ZLINES,
        'Z-order for river lines.'
    ),

    # Subplots settings
    'subplots.align': (
        False,
        'Whether to align axis labels during draw. See `aligning labels '
        '<https://matplotlib.org/stable/gallery/subplots_axes_and_figures/align_labels_demo.html>`__.'  # noqa: E501
    ),
    'subplots.innerpad': (
        '1em',
        'Padding between adjacent subplots. ' + _addendum_units
    ),
    'subplots.outerpad': (
        '0.5em',
        'Padding around figure edge. ' + _addendum_units
    ),
    'subplots.panelpad': (
        '0.5em',
        'Padding between subplots and panels, and between stacked panels. '
        + _addendum_units
    ),
    'subplots.panelwidth': (
        '4em',
        'Width of side panels. ' + _addendum_units
    ),
    'subplots.refwidth': (
        '20em',  # about 3 inches wide
        'Default width of the reference subplot. ' + _addendum_units
    ),
    'subplots.share': (
        3,
        'The axis sharing level, one of ``0``, ``1``, ``2``, or ``3``. '
        'See :ref:`the user guide <ug_share>` for details.'
    ),
    'subplots.span': (
        True,
        'Boolean, toggles spanning axis labels. See `~proplot.ui.subplots` for details.'
    ),
    'subplots.tight': (
        True,
        'Boolean, indicates whether to auto-adjust figure bounds and subplot spacings.'
    ),

    # Super title settings
    'suptitle.color': (
        'black',
        'Figure title color.'
    ),
    'suptitle.size': (
        TITLESIZE,
        'Figure title font size.'
    ),
    'suptitle.weight': (
        'bold',
        'Figure title font weight.'
    ),
    'suptitle.pad': (
        '0.5em',
        'Padding between axes content and the figure super title. '
        + _addendum_units
    ),

    # Text settings
    'text.labelsize': (
        LABELSIZE,
        'Meta setting that changes the label-like sizes '
        '``tick.labelsize``, ``axes.labelsize``, ``legend.fontsize``, '
        "and ``grid.labelsize``. Default is ``'medium'``, i.e. "
        'the value of :rcraw:`font.size`' + _addendum_fonts
    ),
    'text.titlesize': (
        TITLESIZE,
        'Meta setting that changes the title-like sizes '
        '``abc.size``, ``title.size``, ``suptitle.size``, '
        'and row and column label sizes like ``leftlabel.size``. '
        "Default is ``'med-large'``, i.e. 1.1 times :rcraw:`font.size`"
        + _addendum_fonts
    ),

    # Tick settings
    'tick.color': (
        COLOR,
        'Major and minor tick color.'
    ),
    'tick.dir': (
        TICKDIR,
        'Major and minor tick direction. Must be one of '
        "``'out'``, ``'in'``, or ``'inout'``."
    ),
    'tick.labelcolor': (
        COLOR,
        'Axis tick label color. Mirrors the *axis* label '
        ':rcraw:`axes.labelcolor` setting.'
    ),
    'tick.labelsize': (
        LABELSIZE,
        'Axis tick label font size. Mirrors the *axis* label '
        ':rcraw:`axes.labelsize` setting.'
    ),
    'tick.labelweight': (
        'normal',
        'Axis tick label font weight. Mirrors the *axis* label '
        ':rcraw:`axes.labelweight` setting.'
    ),
    'tick.len': (
        TICKLEN,
        'Length of major ticks in points.'
    ),
    'tick.lenratio': (
        TICKLENRATIO,
        'Ratio of minor tickline length to major tickline length.'
    ),
    'tick.minor': (
        TICKMINOR,
        'Boolean, toggles minor ticks on and off.',
    ),
    'tick.pad': (
        TICKPAD,
        'Padding between ticks and tick labels. ' + _addendum_units
    ),
    'tick.ratio': (
        TICKRATIO,
        'Ratio of minor tickline width to major tickline width.'
    ),

    # Title settings
    'title.above': (
        True,
        'Boolean, indicates whether to move outer titles and a-b-c labels above '
        'panels, colorbars, or legends that are above the axes.'
    ),
    'title.pad': (
        TITLEPAD,
        'Padding between the axes edge and the inner and outer titles and '
        'a-b-c labels. Alias for :rcraw:`axes.titlepad`. ' + _addendum_units
    ),
    'title.border': (
        True,
        'Boolean, indicates whether to draw a white border around titles '
        'when :rcraw:`title.loc` is inside the axes.'
    ),
    'title.borderwidth': (
        1.5,
        'Width of the border around titles.'
    ),
    'title.bbox': (
        False,
        'Boolean, whether to draw semi-transparent bounding boxes around titles '
        'when :rcraw:`title.loc` is inside the axes.'
    ),
    'title.bboxcolor': (
        'w',
        'Axes title bounding box color.'
    ),
    'title.bboxstyle': (
        'square',
        'Axes title bounding box style.'
    ),
    'title.bboxalpha': (
        0.5,
        'Axes title bounding box opacity.'
    ),
    'title.bboxpad': (
        None,
        'Padding for the title bounding box. By default this is scaled '
        'to make the box flush against the axes edge. ' + _addendum_units
    ),
    'title.color': (
        'black',
        'Axes title color.'
    ),
    'title.loc': (
        'c',
        'Title position. For options, see the :ref:`title location '
        'table <title_table>`.'
    ),
    'title.size': (
        TITLESIZE,
        'Axes title font size.'
    ),
    'title.weight': (
        'normal',
        'Axes title font weight.'
    ),

    # Top subplot label settings
    'toplabel.color': (
        'black',
        'Font color for column labels on the top of the figure.'
    ),
    'toplabel.size': (
        TITLESIZE,
        'Font size for column labels on the top of the figure.'
    ),
    'toplabel.weight': (
        'bold',
        'Font weight for column labels on the top of the figure.'
    ),
    'toplabel.pad': (
        '0.3em',
        'Padding between axes content and column labels on the top of the figure. '
        + _addendum_units
    ),
}

# Child settings -- changing the parent changes all the children, but changing
# any of the children does not change the parent.
# NOTE: Do not link title.color to axes.titlecolor because the latter
# can have value 'auto' which is not handled in format() right now,
# and because setting was only introduced in version 3.2.
_rc_children = {
    'color': (  # change the 'color' of an axes
        'axes.edgecolor', 'axes.labelcolor',
        'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'
    ),
    'text.labelsize': (  # the 'small' fonts
        'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'grid.labelsize'
    ),
    'text.titlesize': (  # the 'large' fonts
        'abc.size', 'figure.titlesize',
        'axes.titlesize', 'suptitle.size', 'title.size',
        'leftlabel.size', 'toplabel.size',
        'rightlabel.size', 'bottomlabel.size'
    ),
    'linewidth': (
        # NOTE: rc_configurator will adjust [xy]tick.minor.width accordingly
        # NOTE: do not add grid.linewidth to this because common use case is
        # making the edge a bit thicker to *highlight* a subplot, and generally
        # we do not want that to affect gridlines.
        'axes.linewidth', 'xtick.major.width', 'ytick.major.width'
    ),
    'margin': (
        'axes.xmargin', 'axes.ymargin'
    ),
    'tick.color': (
        'xtick.color', 'ytick.color',
    ),
    'tick.dir': (
        'xtick.direction', 'ytick.direction'
    ),
    'tick.len': (
        'xtick.major.size', 'ytick.major.size'
    ),
    'tick.labelsize': (
        'xtick.labelsize', 'ytick.labelsize',
    ),
    'tick.pad': (
        'xtick.major.pad', 'xtick.minor.pad',
        'ytick.major.pad', 'ytick.minor.pad'
    ),
    'grid.color': (
        'gridminor.color', 'grid.labelcolor',
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
    'formatter.limits': (
        'axes.formatter.limits',
    ),
    'formatter.use_locale': (
        'axes.formatter.use_locale',
    ),
    'formatter.use_mathtext': (
        'axes.formatter.use_mathtext',
    ),
    'formatter.min_exponent': (
        'axes.formatter.min_exponent',
    ),
    'formatter.use_offset': (
        'axes.formatter.useoffset',
    ),
    'formatter.offset_threshold': (
        'axes.formatter.offset_threshold',
    ),
}

# Symmetric aliases. Changing one setting changes the other
_rc_aliases = {
    'alpha': 'axes.alpha',
    'axes.titlesize': 'title.size',  # NOTE: translate "auto" to color
    'facecolor': 'axes.facecolor',
    'font.name': 'font.family',
    'grid.below': 'axes.axisbelow',
    'title.pad': 'axes.titlepad',
}
for key, value in _rc_aliases.items():
    _rc_children[key] = (value,)
    _rc_children[value] = (key,)

# Various helper dicts
# NOTE: Make sure to add deprecated rc settings to nodots.
_rc_proplot_default = {
    key: value for key, (value, *_) in _rc_proplot.items()
}

_rc_nodots = {
    name.replace('.', ''): name
    for dict_ in (
        _rc_proplot_default, _rc_matplotlib_default_full, _rc_removed, _rc_renamed
    )
    for name in dict_.keys()
}

_rc_categories = {
    '.'.join(name.split('.')[:i + 1])
    for dict_ in (_rc_proplot_default, _rc_matplotlib_default_full)
    for name in dict_
    for i in range(len(name.split('.')) - 1)
}


def _get_default_param(key):
    """
    Get the default parameter from one of three places. This is
    used for the :rc: role when compiling docs.
    """
    sentinel = object()
    for dict_ in (
        _rc_proplot_default,
        _rc_matplotlib_default,
        _rc_matplotlib_default_full,
    ):
        value = dict_.get(key, sentinel)
        if value is not sentinel:
            return value
    raise KeyError(f'Invalid key {key!r}.')


def _gen_yaml_table(rcdict, comment=True, description=True):
    """
    Return the settings as a nicely tabulated YAML-style table.
    """
    NoneType = type(None)
    prefix = '# ' if comment else ''
    data = []
    for key, pair in rcdict.items():
        # Optionally append description
        if description and isinstance(pair, tuple):  # add commented out description
            value = pair[0]
            descrip = '# ' + pair[1]
        else:
            descrip = ''
            if isinstance(pair, tuple):
                value = pair[0]
            else:
                value = pair

        # Translate object to string
        if isinstance(value, cycler.Cycler):  # special case!
            value = repr(value)
        elif isinstance(value, (str, numbers.Number, NoneType)):
            value = str(value)
        elif isinstance(value, (list, tuple, np.ndarray)) and all(
            isinstance(val, (str, numbers.Number)) for val in value
        ):
            value = ', '.join(str(val) for val in value)
        else:
            warnings._warn_proplot(
                f'Failed to write rc setting {key} = {value!r}. Must be string, '
                'number, or list or tuple thereof, or None or a cycler.'
            )
            continue
        if value[:1] == '#':  # e.g. HEX string
            value = repr(value)
        data.append((key, value, descrip))

    # Generate string
    string = ''
    keylen = len(max(rcdict, key=len))
    vallen = len(max((tup[1] for tup in data), key=len))
    for key, value, descrip in data:
        space1 = ' ' * (keylen - len(key) + 1)
        space2 = ' ' * (vallen - len(value) + 2) if descrip else ''
        string += f'{prefix}{key}:{space1}{value}{space2}{descrip}\n'

    return string.strip()


def _gen_rst_table():
    """
    Return the settings in an RST-style table.
    """
    # Initial stuff
    colspace = 2  # spaces between each column
    descrips = tuple(descrip for _, descrip in _rc_proplot.values())
    keylen = len(max((*_rc_proplot, 'Key'), key=len)) + 4  # for literal backticks
    vallen = len(max((*descrips, 'Description'), key=len))
    divider = '=' * keylen + ' ' * colspace + '=' * vallen + '\n'
    header = 'Key' + ' ' * (keylen - 3 + colspace) + 'Description\n'

    # Build table
    string = divider + header + divider
    for key, (_, descrip) in _rc_proplot.items():
        spaces = ' ' * (keylen - (len(key) + 4) + colspace)
        string += f'``{key}``{spaces}{descrip}\n'

    string = string + divider
    return '.. rst-class:: proplot-rctable\n\n' + string.strip()
