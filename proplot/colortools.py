#!/usr/bin/env python3
"""
Registers colormaps, color cycles, and color string names with
`register_cmaps`, `register_cycles`, and `register_colors`.
Defines tools for creating new colormaps and color cycles, i.e. `colormap`
and `colors`. Defines helpful new `~matplotlib.colors.Normalize` and
`~matplotlib.colors.Colormap` classes.

For a visual reference, see the :ref:`Table of colormaps`,
:ref:`Table of color cycles`, and the :ref:`Table of colors`.

Perceptually uniform colormaps
------------------------------

ProPlot's custom colormaps are instances of the new
`PerceptuallyUniformColormap` class. These classes employ *linear transitions*
between channel values in any of three possible "perceptually uniform",
HSV-like colorspaces.  These colorspaces can be described as follows (see also
`these visualizations <http://www.hsluv.org/comparison/>`_):

* **HCL**: A purely perceptually uniform colorspace, where colors are
  broken down into “hue” (color, range 0-360), “chroma”
  (colorfulness, range 0-100), and “luminance” (brightness, range 0-100).
* **HPLuv**: As with HCL, but 100 chroma is scaled to be the *minimum maximum
  chroma* across all hues for a given luminance, and is hence more
  appropriate for multi-hue colormaps.
* **HSLuv**: As with HCL, but 100 chroma is scaled to be the *maximum
  possible chroma* for a given hue and luminance. This is more appropriate for
  single-hue colormaps, because crossing hues in
  this space make it more likely that bands of higher absolute chroma are
  crossed.

The HCL space is the only "purely" perceptually uniform colorspace. But
during a linear transition between two values, we may cross over "impossible"
colors (i.e. colors with RGB channels >1).

The HSLuv and HPLuv colorspaces
were developed to resolve this issue by (respectively) scaling and clipping
high-chroma colors across different hues and luminances.

From other projects
-------------------

I’ve removed some outdated “miscellaneous” colormaps that are packaged
by default (see `this reference
<https://matplotlib.org/examples/color/colormaps_reference.html>`_)
and added new "perceptually uniform" colormaps from the following projects:

* The `cmOcean project <https://matplotlib.org/cmocean/>`_
* The `SciVisColor project <https://sciviscolor.org/home/colormaps/>`_
* `Kenneth Moreland's colormaps <http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/index.html>`_
* `Fabio Crameri's colormaps <http://www.fabiocrameri.ch/colourmaps.php>`_
* Colormaps commissioned by `Statistik Stadt Zürich
  <http://soliton.vm.bytemark.co.uk/pub/cpt-city/ssz/index.html>`_
* Peter Koveski's `CET colormaps <https://peterkovesi.com/projects/colourmaps/>`_

Several of these were found thanks to `Riley X. Bradey
<https://github.com/bradyrx>`_. Others were found using the `cpt-city
<http://soliton.vm.bytemark.co.uk/pub/cpt-city/>`_ archive of color
gradients.

Note that matplotlib comes packaged with every `ColorBrewer2.0
<http://colorbrewer2.org/>`__ colormap, which are also certainly
"perceptually uniform".


Flexible colormap arguments
---------------------------

All of the `~matplotlib.axes.Axes` methods listed in
`~proplot.axes.cmap_methods` and
`~proplot.axes.cycle_methods` have been wrapped by ProPlot. For the
latter methods, ProPlot adds a brand new keyword arg called ``cycle``, used for
changing the axes property cycler on-the-fly.

The ``cmap`` and ``cycle`` arguments
are all passed through the magical `Colormap` function.
`Colormap` is incredibly powerful -- it can make colormaps
on-the-fly, look up existing maps, and merge them. As such, any of the following
are now valid ``cmap`` and ``cycle`` arguments:

1. Registered colormap names. For example, ``'Blues'`` or ``'Sunset'``.
   See :ref:`Table of colormaps` for a visalization of the registered maps.
   Cycles are constructed automatically be sampling colors from these colormaps;
   use e.g. ``('Blues', 5)`` to specify the number of colors in the cycle.
2. Gradations of a single hue. For example, ``'maroon'`` creates colormap spanning
   white to maroon, and ``'maroon90'`` spans from a 90% luminance pale red
   color to maroon (see `~proplot.colors.Colormap`).
   Cycles are again constructed automatically; use e.g. ``('maroon', 10)``
   to specify the number of colors in the cycle.
3. Registered color cycle names. For example, ``'Cycle1'``. See
   :ref:`Table of color cycles`
   for a visualization of the registered cycles.
4. Lists of colors. For example, ``['red', 'blue', 'green']``.
   This makes a `~matplotlib.colors.ListedColormap` map which can be trivially
   used as a "color cycle".
5. Dictionary containing the keys ``'h'``, ``'s'``, and ``'l'``. This builds
   a `PerceptuallyUniformColormap` using the
   `~PerceptuallyUniformColormap.from_hsl` constructor.
6. **List of any of the above five arguments, to merge the resulting
   colormaps.** For
   example, use ``['viridis', 'navy']`` to merge the virids map with a colormap
   spanning white to navy.

Note when assigning to the ``proplot.rc.cycle`` global setting (see
`~proplot.rcmod`), the argument is also interpreted as above. For example,
``proplot.rc.cycle = ('blue', 10)`` will construct a color cycle with 10 colors
ranging from white to blue.
"""
#------------------------------------------------------------------------------#
# Interesting cpt-city colormaps that did not use:
# * Considered Jim Mossman maps, but not very uniform.
# * Erik Jeschke grayscale ones are also neat, but probably not much
#   scientific use.
# * Rafi 'sky' themes were pretty, but ultimately not useful.
# * Crumblingwalls also interesting, but too many/some are weird.
# * NCL gradients mostly ugly, forget them.
# * Piecrust design has interesting 'nature' colormaps, but they are
#   not practical. Just added a couple astro ones (aurora, space, star).
# * Elvensword is cool, but most are pretty banded.
# Geographic ones not used:
# * Christian Heine geologic time maps are interesting, but again not
#   uniform and not useful.
# * IBCA could have been good, but bathymetry has ugly jumps.
# * GMT maps also interesting, but non uniform.
# * Victor Huérfano Caribbean map almost useful, but has banding.
# * Christopher Wesson martian topo map sort of cool, but too pale.
# * Thomas Deweez has more topo colors, but kind of ugly.
# Geographic ones to be adapted:
# * ESRI seems to have best geographic maps.
# * Wiki schemes are also pretty good, but not great.
#------------------------------------------------------------------------------#
# Notes on 'channel-wise alpha':
# * Colormaps generated from HCL space (and cmOcean ones) are indeed perfectly
#   perceptually uniform, but this still looks bad sometimes -- usually we
#   *want* to focus on the *extremes*, so want to weight colors more heavily
#   on the brighters/whiter part of the map! That's what the ColdHot map does,
#   it's what most of the ColorBrewer maps do, and it's what ColorWizard does.
# * By default extremes are stored at end of *lookup table*, not as
#   separate RGBA values (so look under cmap._lut, indexes cmap._i_over and
#   cmap._i_under). You can verify that your cmap is using most extreme values
#   by comparing high-resolution one to low-resolution one.
#------------------------------------------------------------------------------#
# Potential bottleneck, loading all this stuff?
# NO. Try using @timer on register functions, turns out worst is colormap
# one at 0.1 seconds. Just happens to be a big package, takes a bit to compile
# to bytecode (done every time module changed) then import.
#------------------------------------------------------------------------------#
# Here's some useful info on colorspaces
# https://en.wikipedia.org/wiki/HSL_and_HSV
# http://www.hclwizard.org/color-scheme/
# http://www.hsluv.org/comparison/ compares lch, hsluv (scaled lch), and hpluv (truncated lch)
# Info on the CIE conventions
# https://en.wikipedia.org/wiki/CIE_1931_color_space
# https://en.wikipedia.org/wiki/CIELUV
# https://en.wikipedia.org/wiki/CIELAB_color_space
# And some useful tools for creating colormaps and cycles
# https://nrlmry.navy.mil/TC.html
# http://help.mail.colostate.edu/tt_o365_imap.aspx
# http://schumacher.atmos.colostate.edu/resources/archivewx.php
# https://coolors.co/
# http://tristen.ca/hcl-picker/#/hlc/12/0.99/C6F67D/0B2026
# http://gka.github.io/palettes/#diverging|c0=darkred,deeppink,lightyellow|c1=lightyellow,lightgreen,teal|steps=13|bez0=1|bez1=1|coL0=1|coL1=1
# https://flowingdata.com/tag/color/
# http://tools.medialab.sciences-po.fr/iwanthue/index.php
# https://learntocodewith.me/posts/color-palette-tools/
#------------------------------------------------------------------------------
import os
import re
import json
import glob
from lxml import etree
from numbers import Number
import cycler
import warnings
import numpy as np
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from matplotlib import rcParams # cannot import rcmod because rcmod import this!
from . import utils, colormath
from .utils import _default, _counter, ic
_data_user = os.path.join(os.path.expanduser('~'), '.proplot')
_data_cmaps = os.path.join(os.path.dirname(__file__), 'cmaps') # or parent, but that makes pip install distribution hard
_data_colors = os.path.join(os.path.dirname(__file__), 'colors') # or parent, but that makes pip install distribution hard

# Define some new palettes
# Note the default listed colormaps
_cycles_loaded = {}
_cycles_preset = {
    # Default matplotlib v2
    'default':      ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    # Copied from stylesheets; stylesheets just add color themese from every possible tool, not already present as a colormap
    '538':          ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'],
    'ggplot':       ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
    # The default seaborn ones (excluded deep/muted/bright because thought they were unappealing)
    'ColorBlind':   ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9'],
    'ColorBlind10': ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"], # versions with more colors
    # From the website
    'FlatUI':       ["#3498db", "#e74c3c", "#95a5a6", "#34495e", "#2ecc71", "#9b59b6"],
    # Created with online tools; add to this
    # See: http://tools.medialab.sciences-po.fr/iwanthue/index.php
    'Cinematic':    [(51,92,103), (158,42,43), (255,243,176), (224,159,62), (84,11,14)],
    'Cool':         ["#6C464F", "#9E768F", "#9FA4C4", "#B3CDD1", "#C7F0BD"],
    'Sugar':        ["#007EA7", "#B4654A", "#80CED7", "#B3CDD1", "#003249"],
    'Vibrant':      ["#007EA7", "#D81159", "#B3CDD1", "#FFBC42", "#0496FF"],
    'Office':       ["#252323", "#70798C", "#DAD2BC", "#F5F1ED", "#A99985"],
    'Industrial':   ["#38302E", "#6F6866", "#788585", "#BABF95", "#CCDAD1"],
    'Tropical':     ["#0D3B66", "#F95738", "#F4D35E", "#FAF0CA", "#EE964B"],
    'Intersection': ["#2B4162", "#FA9F42", "#E0E0E2", "#A21817", "#0B6E4F"],
    'Field':        ["#23395B", "#D81E5B", "#FFFD98", "#B9E3C6", "#59C9A5"],
    }

# Color stuff
# Keep major color names, and combinations of those names
_distinct_colors_threshold = 0.09 # bigger number equals fewer colors
_distinct_colors_space = 'hcl' # register colors distinct in this space?
_exceptions_names = (
    'sky blue', 'eggshell', 'sea blue', 'coral', 'tomato red', 'brick red', 'crimson',
    'red orange', 'yellow orange', 'yellow green', 'blue green',
    'blue violet', 'red violet',
    )
_bad_names = '(' + '|'.join(( # filter these out; let's try to be professional here...
    'shit', 'poo', 'pee', 'piss', 'puke', 'vomit', 'snot', 'booger',
    )) + ')'
_sanitize_names = ( # replace regex (first entry) with second entry
    ('/', ' '), ("'s", ''), ('grey', 'gray'),
    ('pinky', 'pink'), ('greeny', 'green'),
    ('bluey',  'blue'),
    ('robin egg', 'robins egg'),
    ('egg blue', 'egg'),
    (r'reddish', 'red'),
    (r'purplish', 'purple'),
    (r'bluish',  'blue'),
    (r'ish\b', ''),
    ('bluegray', 'blue gray'),
    ('grayblue', 'gray blue'),
    ('lightblue', 'light blue')
    )
_space_aliases = {
    'rgb':   'rgb',
    'hsv':   'hsv',
    'hpl':   'hpl',
    'hpluv': 'hpl',
    'hsl':   'hsl',
    'hsluv': 'hsl',
    'hcl':   'hcl',
    'lch':   'hcl',
    }
_channel_idxs = {
    'h': 0, 'hue': 0,
    's': 1, 'saturation': 1,
    'c': 1, 'chroma': 1,
    'l': 2, 'luminance': 2,
    'a': 3, 'alpha': 3,
    }

# Names of builtin colormaps
# NOTE: Has support for 'x' coordinates in first column.
# NOTE: For 'alpha' column, must use a .rgba filename
# TODO: Better way to save colormap files.
_cmap_categories = {
    # Your custom registered maps; this is a placeholder, meant to put these
    # maps at the top of the colormap table
    'User': [
        ],

    # We keep these ones
    'Matplotlib Originals': [
        'viridis', 'plasma', 'inferno', 'magma', 'twilight', 'twilight_shifted',
        ],

    # Assorted origin, but these belong together
    'Grayscale': [
        'Grays',
        'Mono',
        'GrayCycle',
        'GrayCycle_shifted',
        ],

    # CET isoluminant maps
    # See: https://peterkovesi.com/projects/colourmaps/
    # All the others have better options
    'Isoluminant': [
        'Iso1', 'Iso2', 'Iso3', 
        'Phase', 'Phase_shifted', # these actually from cmocean
        ],

    # Included ColorBrewer
    'ColorBrewer2.0 Sequential': [
        'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'PuBu', 'PuBuGn', 'BuGn', 'GnBu', 'YlGnBu', 'YlGn'
        ],

    # Added diverging versions
    # See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/jjg/polarity/index.html
    # Other JJ Green maps weren't that great
    # TODO: Add 'serated' maps? See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/jjg/serrate/index.html
    # TODO: Add tool for cutting center out of ***any*** colormap by ending
    # with the _cut suffix or something?
    'ColorBrewer2.0 Diverging': [
        'Spectral', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGY',
        'RdBu', 'RdYlBu', 'RdYlGn',
        ],

    # Custom maps
    'ProPlot Sequential': [
         'Glacial',
        'Bog', 'Verdant',
        'Turquoise',
        'Sunrise', 'Sunset', 'Fire',
        'Golden'
        ],
        # 'Vibrant'], # empty at first, fill automatically
    'ProPlot Diverging': [
        'IceFire', 'NegPos', 'BlueRed', 'PurplePink', 'DryWet', 'AltDryWet', 'LandSea'
        ],

    # Other
    # BlackBody2 is actually from Los Alamos, and a couple are from Kenneth's
    # website, but organization is better this way.
    'Misc Diverging': [
        'bwr',
        'CoolWarm',
        'SoftCoolWarm',
        'MutedCoolWarm',
        'ColdHot',
        'Temp', # from ???
        # 'BlackBody1', 'BlackBody2', 'BlackBody3', # 3rd one is actually sky theme from rafi
        # 'Star',
        # 'JMN', # James map; ugly, so deleted
        # 'CubeHelix', 'SatCubeHelix',
        # 'cividis',
        # 'Aurora', 'Space', # from PIEcrust; not uniform, so deleted
        # 'TemperatureJJG', # from JJG; ugly, so deleted
        # 'Kindlmann', 'ExtendedKindlmann',
        # 'Seismic', # note this one originally had hard boundaries/no interpolation
        # 'MutedBio', 'DarkBio', # from: ???, maybe SciVisColor
        ],

    # Statistik
    # 'Statistik Stadt Zürich': [
    'Zürich Muted': [
        'MutedBlue', 'MutedRed', 'MutedDry', 'MutedWet',
        'MutedBuRd', 'MutedBuRd_cut', 'MutedDryWet', 'MutedDryWet_cut',
        ],

    # cmOcean
    'cmOcean Sequential': [
        'Oxy', 'Thermal', 'Dense', 'Ice', 'Haline',
        'Deep', 'Algae', 'Tempo', 'Speed', 'Turbid', 'Solar', 'Matter',
        'Amp',
        ],
    'cmOcean Diverging': [
        'Balance', 'Curl', 'Delta'
        ],

    # SciVisColor
    # Culled these because some were ugly
    # Actually nevermind... point of these is to *combine* them, make
    # stacked colormaps that highlight different things.
    'SciVisColor Blues': [
        'Blue0', 'Blue1', 'Blue2', 'Blue3', 'Blue4', 'Blue5', 'Blue6', 'Blue7', 'Blue8', 'Blue9', 'Blue10', 'Blue11',
        ],
    'SciVisColor Greens': [
        'Green1', 'Green2', 'Green3', 'Green4', 'Green5', 'Green6', 'Green7', 'Green8',
        ],
    'SciVisColor Oranges': [
        'Orange1', 'Orange2', 'Orange3', 'Orange4', 'Orange5', 'Orange6', 'Orange7', 'Orange8',
        ],
    'SciVisColor Browns': [
        'Brown1', 'Brown2', 'Brown3', 'Brown4', 'Brown5', 'Brown6', 'Brown7', 'Brown8', 'Brown9',
        ],
    'SciVisColor Reds/Purples': [
        'RedPurple1', 'RedPurple2', 'RedPurple3', 'RedPurple4', 'RedPurple5', 'RedPurple6', 'RedPurple7', 'RedPurple8',
        ],

    # FabioCrameri
    # See: http://www.fabiocrameri.ch/colourmaps.php
    'Fabio Crameri Sequential': [
        'Acton', 'Buda', 'Lajolla',
        'Imola', 'Bamako', 'Nuuk', 'Davos', 'Oslo', 'Devon', 'Tokyo', 'Hawaii', 'Batlow',
        'Turku', 'Bilbao', 'Lapaz',
        ],
    'Fabio Crameri Diverging': [
        'Roma', 'Broc', 'Cork',  'Vik', 'Oleron', 'Lisbon', 'Tofino', 'Berlin',
        ],

    # Gross. These ones will be deleted.
    'Alt Sequential': [
        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'multi', 'cividis',
        'afmhot', 'gist_heat', 'copper'
        ],
    'Alt Rainbow': [
        'multi', 'cividis'
        ],
    'Alt Diverging': [
        'coolwarm', 'bwr', 'seismic'
        ],
    'Miscellaneous Orig': [
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'brg', 'hsv', 'hot', 'rainbow',
        'gist_rainbow', 'jet', 'nipy_spectral', 'gist_ncar', 'cubehelix',
        ],

    # Kenneth Moreland
    # See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/index.html
    # Soft coolwarm from: https://www.kennethmoreland.com/color-advice/
    # 'Kenneth Moreland': [
    #     'CoolWarm', 'MutedCoolWarm', 'SoftCoolWarm',
    #     'BlueTan', 'PurpleOrange', 'CyanMauve', 'BlueYellow', 'GreenRed',
    #     ],
    # 'Kenneth Moreland Sequential': [
    #     'BlackBody', 'Kindlmann', 'ExtendedKindlmann',
    #     ],

    # Los Alamos
    # See: https://datascience.lanl.gov/colormaps.html
    # Most of these have analogues in SciVisColor, previously added the few
    # unique ones to Miscellaneous category
    # 'Los Alamos Sequential': [
    #     'MutedRainbow', 'DarkRainbow', 'MutedBlue', 'DeepBlue', 'BrightBlue', 'BrightGreen', 'WarmGray',
    #     ],
    # 'Los Alamos Diverging': [
    #     'MutedBlueGreen', 'DeepBlueGreen', 'DeepBlueGreenAsym', 'DeepColdHot', 'DeepColdHotAsym', 'ExtendedCoolWarm'
    #     ],

    # Removed the "combo" maps (and ugly diverging ones) because these can
    # be built in proplot with the Colormap tool!
    # 'SciVisColor Diverging': [
    #     'Div1', 'Div2', 'Div3', 'Div4', 'Div5'
    #     ],
    # 'SciVisColor 3 Waves': [
    #     '3Wave1', '3Wave2', '3Wave3', '3Wave4', '3Wave5', '3Wave6', '3Wave7'
    #     ],
    # 'SciVisColor 4 Waves': [
    #     '4Wave1', '4Wave2', '4Wave3', '4Wave4', '4Wave5', '4Wave6', '4Wave7'
    #     ],
    # 'SciVisColor 5 Waves': [
    #     '5Wave1', '5Wave2', '5Wave3', '5Wave4', '5Wave5', '5Wave6'
    #     ],
    # 'SciVisColor Waves': [
    #     '3Wave1', '3Wave2', '3Wave3',
    #     '4Wave1', '4Wave2', '4Wave3',
    #     '5Wave1', '5Wave2', '5Wave3',
    #     ],
    # 'SciVisColor Inserts': [
    #     'Insert1', 'Insert2', 'Insert3', 'Insert4', 'Insert5', 'Insert6', 'Insert7', 'Insert8', 'Insert9', 'Insert10'
    #     ],
    # 'SciVisColor Thick Inserts': [
    #     'ThickInsert1', 'ThickInsert2', 'ThickInsert3', 'ThickInsert4', 'ThickInsert5'
    #     ],
    # 'SciVisColor Highlight': [
    #     'Highlight1', 'Highlight2', 'Highlight3', 'Highlight4', 'Highlight5',
    #     ],

    # Most of these were ugly, deleted them
    # 'SciVisColor Outlier': [
    #     'DivOutlier1', 'DivOutlier2', 'DivOutlier3', 'DivOutlier4',
    #     'Outlier1', 'Outlier2', 'Outlier3', 'Outlier4'
    #     ],

    # Duncan Agnew
    # See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/dca/index.html
    # These are 1.0.5 through 1.4.0
    # 'Duncan Agnew': [
    #     'Alarm1', 'Alarm2', 'Alarm3', 'Alarm4', 'Alarm5', 'Alarm6', 'Alarm7'
    #     ],

    # Elevation and bathymetry
    # 'Geographic': [
    #     'Bath1', # from XKCD; see: http://soliton.vm.bytemark.co.uk/pub/cpt-city/xkcd/tn/xkcd-bath.png.index.html
    #     'Bath2', # from Tom Patterson; see: http://soliton.vm.bytemark.co.uk/pub/cpt-city/tp/index.html
    #     'Bath3', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/ibcso/tn/ibcso-bath.png.index.html
    #     'Bath4', # ^^ same
    #     'Geography4-1', # mostly ocean
    #     'Geography5-4', # range must be -4000 to 5000
    #     'Geography1', # from ???
    #     'Geography2', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/ngdc/tn/ETOPO1.png.index.html
    #     'Geography3', # from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/mby/tn/mby.png.index.html
    #     ],
    }
# Categories to ignore/*delete* from dictionary because they suck donkey balls
_cmap_categories_delete = ['Alt Diverging', 'Alt Sequential', 'Alt Rainbow', 'Miscellaneous Orig']

# Slice indices that split up segments of names
# WARNING: Must add to this list manually! Not worth trying to generalize.
# List of string cmap names, and the indices where they can be broken into parts
_cmap_parts = {
    # Diverging colorbrewer
    'piyg':         (None, 2, None),
    'prgn':         (None, 1, 2, None), # purple red green
    'brbg':         (None, 2, 3, None), # brown blue green
    'puor':         (None, 2, None),
    'rdgy':         (None, 2, None),
    'rdbu':         (None, 2, None),
    'rdylbu':       (None, 2, 4, None),
    'rdylgn':       (None, 2, 4, None),
    # Other diverging
    'coldhot':      (None, 4, None),
    'bwr':          (None, 1, 2, None),
    'icefire':      (None, 3, None),
    'negpos':       (None, 3, None),
    'bluered':      (None, 4, None),
    'purplepink':   (None, 4, None),
    'drywet':       (None, 3, None),
    'drierwetter':  (None, 5, None),
    'landsea':      (None, 4, None),
    }
# Tuple pairs of mirror image cmap names
_cmap_mirrors = [
    (name, ''.join(reversed([name[slice(*idxs[i:i+2])] for i in range(len(idxs)-1)])),)
    for name,idxs in _cmap_parts.items()
    ]

#------------------------------------------------------------------------------#
# Special classes
#------------------------------------------------------------------------------#
# Class for flexible color names
# WARNING: Matplotlib 'color' arguments are passed to to_rgba, which tries
# to read directly from cache and if that fails, tries to sanitize input.
# The sanitization raises error when encounters (colormap, idx) tuple. So
# we need to override the *cache* instead of color dictionary itself!
# WARNING: Builtin to_rgb tries to get cached colors as dict[name, alpha],
# resulting in key as (colorname, alpha) or ((R,G,B), alpha) tuple. Impossible
# to differentiate this from (cmapname, index) usage! Must do try except lookup
# into colormap dictionary every time. Don't want to do this for actual
# color dict for sake of speed, so we only wrap *cache* lookup. Also we try
# to avoid cmap lookup attempt whenever possible.
class ColorDictSpecial(dict):
    """Special dictionary that lets user draw single color tuples from
    arbitrary colormaps or color cycles."""
    def __getitem__(self, key):
        """
        Either samples the color from a colormap or color cycle,
        or calls the parent getitem to look up the color name.

        For a **smooth colormap**, usage is e.g.
        ``color=('Blues', 0.8)`` -- the number should be between 0 and 1, and
        indicates where to draw the color from the smooth colormap. For a
        "listed" colormap, i.e. a **color cycle**, usage is e.g.
        ``color=('colorblind', 2)``. The number indicates the index in the
        list of discrete colors.

        These examples work with any matplotlib command that accepts
        a ``color`` keyword arg.
        """
        # Pull out alpha
        # WARNING: Possibly fragile? Does this hidden behavior ever change?
        # NOTE: This override doubles startup time 0.0001s to 0.0002s, probably ok.
        if np.iterable(key) and len(key)==2:
            key, alpha = key
        if np.iterable(key) and len(key)==2 and \
            isinstance(key[1], Number) and isinstance(key[0], str): # i.e. is not None; this is *very common*, so avoids lots of unnecessary lookups!
            try:
                cmap = mcm.cmap_d[key[0]]
            except (TypeError, KeyError):
                pass
            else:
                if isinstance(cmap, mcolors.ListedColormap):
                    return tuple(cmap.colors[key[1]]) # draw color from the list of colors, using index
                else:
                    return tuple(cmap(key[1])) # interpolate color from colormap, using key in range 0-1
        return super().__getitem__((key, alpha))
# Wraps the cache
class _ColorMappingOverride(mcolors._ColorMapping):
    def __init__(self, mapping):
        """Wraps the cache."""
        super().__init__(mapping)
        self.cache = ColorDictSpecial({})
# Override default color name dictionary
if not isinstance(mcolors._colors_full_map, _ColorMappingOverride):
    mcolors._colors_full_map = _ColorMappingOverride(mcolors._colors_full_map)

# List of colors with 'name' attribute
class CycleList(list):
    """Simply stores a list of colors, and adds a `name` attribute corresponding
    to the registered name."""
    def __repr__(self):
        """Wraps the string representation."""
        return 'CycleList(' + super().__repr__() + ')'
    def __init__(self, list_, name):
        self.name = name
        super().__init__(list_)

# Flexible colormap identification
class CmapDict(dict):
    """
    Flexible, case-insensitive colormap identification. Replaces the
    `matplotlib.cm.cmap_d` dictionary that stores registered colormaps.

    Behaves like a dictionary, with three new features:

    1. Names are case insensitive: ``'Blues'``, ``'blues'``, and ``'bLuEs'``
       are all valid names for the "Blues" colormap.
    2. "Reversed" colormaps are not stored directly: Requesting e.g.
       ``'Blues_r'`` will just look up ``'Blues'``, then return the result
       of the `~matplotlib.colors.Colormap.reversed` method.
    3. Diverging colormap names can be referenced by their "inverse" name.
       For example, ``'BuRd'`` is equivalent to ``'RdBu_r'``, as are
       ``'BuYlRd'`` and ``'RdYlBu_r'``.
    """
    def __init__(self, kwargs):
        kwargs_filtered = {}
        for key,value in kwargs.items():
            if not isinstance(key, str):
                raise KeyError(f'Invalid key {key}. Must be string.')
            if key[-2:] != '_r': # don't need to store these!
                kwargs_filtered[key.lower()] = value
        super().__init__(kwargs_filtered)

    # Helper functions
    def _sanitize_key(self, key):
        """Sanitizes key name."""
        # Try retrieving
        if not isinstance(key, str):
            raise ValueError(f'Invalid key {key}. Must be string.')
        key = key.lower()
        reverse = False
        if key[-2:] == '_r':
            key = key[:-2]
            reverse = True
        if not super().__contains__(key):
            # Attempt to get 'mirror' key, maybe that's the one
            # stored in colormap dict
            key_mirror = key
            for mirror in _cmap_mirrors:
                try:
                    idx = mirror.index(key)
                    key_mirror = mirror[1 - idx]
                except ValueError:
                    continue
            if super().__contains__(key_mirror):
                reverse = (not reverse)
                key = key_mirror
        # Return 'sanitized' key. Not necessarily in dictionary! Error
        # will be raised further down the line if so.
        if reverse:
            key = key + '_r'
        return key

    def _getitem(self, key):
        """Get value, but skip key sanitization."""
        reverse = False
        if key[-2:] == '_r':
            key = key[:-2]
            reverse = True
        value = super().__getitem__(key) # may raise keyerror
        if reverse:
            try:
                value = value.reversed()
            except AttributeError:
                raise KeyError(f'Dictionary value in {key} must have reversed() method.')
        return value

    # Indexing and 'in' behavior
    def __getitem__(self, key):
        """Sanitizes key, then queries dictionary."""
        # Assume lowercase
        key = self._sanitize_key(key)
        return self._getitem(key)

    def __setitem__(self, key, item):
        """Assigns lowercase."""
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key}. Must be string.')
        return super().__setitem__(key.lower(), item)

    def __contains__(self, item):
        """The 'in' behavior."""
        try:
            self.__getitem__(item)
            return True
        except KeyError:
            return False

    # Other methods
    def get(self, key, *args):
        """Case-insensitive version of `dict.get`."""
        # Get item
        if len(args)>1:
            raise ValueError(f'Accepts only 1-2 arguments (got {len(args)+1}).')
        try:
            if not isinstance(key, str):
                raise KeyError(f'Invalid key {key}. Must be string.')
            return self.__getitem__(key.lower())
        except KeyError as key_error:
            if args:
                return args[0]
            else:
                raise key_error

    def pop(self, key, *args):
        """Case-insensitive version of `dict.pop`."""
        # Pop item
        if len(args)>1:
            raise ValueError(f'Accepts only 1-2 arguments (got {len(args)+1}).')
        try:
            key = self._sanitize_key(key)
            value = self._getitem(key) # could raise error
            del self[key]
        except KeyError as key_error:
            if args:
                return args[0]
            else:
                raise key_error
        return value
# Override default colormap dictionary
if not isinstance(mcm.cmap_d, CmapDict):
    mcm.cmap_d = CmapDict(mcm.cmap_d)

#------------------------------------------------------------------------------#
# Color manipulation functions
#------------------------------------------------------------------------------#
def _get_space(space):
    """Verify requested colorspace is valid."""
    space = _space_aliases.get(space, None)
    if space is None:
        raise ValueError(f'Unknown colorspace "{space}".')
    return space

def _get_channel(color, channel, space='hsl'):
    """Gets hue, saturation, or luminance channel value from registered
    string color name. The color name `color` can optionally be a string
    with the format ``'color+x'`` or ``'color-x'``, where `x` specifies
    the offset from the channel value."""
    # Interpret channel
    channel = _channel_idxs.get(channel, channel)
    if callable(color) or isinstance(color, Number):
        return color
    if channel not in (0,1,2):
        raise ValueError('Channel must be in [0,1,2].')
    # Interpret string or RGB tuple
    offset = 0
    if isinstance(color, str):
        regex = '([-+]\S*)$' # user can optionally offset from color; don't filter to just numbers, want to raise our own error if user messes up
        match = re.search(regex, color)
        if match:
            try:
                offset = float(match.group(0))
            except ValueError:
                raise ValueError(f'Invalid channel identifier "{color}".')
            color = color[:match.start()]
    return offset + to_xyz(to_rgb(color, 'rgb'), space)[channel]

def shade(color, shade=0.5):
    """Changes the "shade" of a color by scaling its luminance channel by `shade`."""
    try:
        color = mcolors.to_rgb(color) # ensure is valid color
    except Exception:
        raise ValueError(f'Invalid RGBA argument {color}. Registered colors are: {", ".join(mcolors._colors_full_map.keys())}.')
    color = [*colormath.rgb_to_hsl(*color)]
    color[2] = max([0, min([color[2]*shade, 100])]) # multiply luminance by this value
    color = [*colormath.hsl_to_rgb(*color)]
    return tuple(color)

def to_rgb(color, space='rgb'):
    """Generalization of matplotlib's `~matplotlib.colors.to_rgb`. Translates
    colors from *any* colorspace to rgb. Also will convert color
    strings to tuple. Inverse of `to_xyz`."""
    # First the RGB input
    # NOTE: Need isinstance here because strings stored in numpy arrays
    # are actually subclasses thereof!
    if isinstance(color, str):
        try:
            color = mcolors.to_rgb(color) # ensure is valid color
        except Exception:
            raise ValueError(f'Invalid RGBA argument {color}. Registered colors are: {", ".join(mcolors._colors_full_map.keys())}.')
    elif space=='rgb':
        color = color[:3] # trim alpha
        if any(c>1 for c in color):
            color = [c/255 for c in color] # scale to within 0-1
    # Next the perceptually uniform versions
    elif space=='hsv':
        color = colormath.hsl_to_rgb(*color)
    elif space=='hpl':
        color = colormath.hpluv_to_rgb(*color)
    elif space=='hsl':
        color = colormath.hsluv_to_rgb(*color)
    elif space=='hcl':
        color = colormath.hcl_to_rgb(*color)
    else:
        raise ValueError('Invalid RGB value.')
    return color

def to_xyz(color, space):
    """Translates from RGB space to colorspace `space`. Inverse of `to_rgb`."""
    # Run tuple conversions
    # NOTE: Don't pass color tuple, because we may want to permit out-of-bounds RGB values to invert conversion
    if isinstance(color, str):
        color = mcolors.to_rgb(color) # convert string
    else:
        color = color[:3]
    if space=='hsv':
        color = colormath.rgb_to_hsl(*color) # rgb_to_hsv would also work
    elif space=='hpl':
        color = colormath.rgb_to_hpluv(*color)
    elif space=='hsl':
        color = colormath.rgb_to_hsluv(*color)
    elif space=='hcl':
        color = colormath.rgb_to_hcl(*color)
    elif space=='rgb':
        pass
    else:
        raise ValueError(f'Invalid colorspace {space}.')
    return color

#------------------------------------------------------------------------------#
# Helper functions
#------------------------------------------------------------------------------#
def _transform_cycle(color):
    """Transforms colors C0, C1, etc. into their corresponding color strings.
    May be necessary trying to change the color cycler."""
    # Optional exit
    if not isinstance(color, str):
        return color
    elif not re.match('^C[0-9]$', color):
        return color
    # Transform color to actual cycle color
    else:
        cycler = rcParams['axes.prop_cycle'].by_key()
        if 'color' not in cycler:
            cycle = ['k']
        else:
            cycle = cycler['color']
        return cycle[int(color[-1])]

def _clip_colors(colors, mask=True, gray=0.2, verbose=False):
    """
    Clips impossible colors rendered in an HSl-to-RGB colorspace conversion.
    Used by `PerceptuallyUniformColormap`. If `mask` is ``True``, impossible
    colors are masked out

    Parameters
    ----------
    colors : list of length-3 tuples
        The RGB colors.
    mask : bool, optional
        Whether to mask out (set to `gray` color) or clip (limit
        range of each channel to 0-1) the out-of-range RGB channels.
    gray : float, optional
        The identical RGB channel values (gray color) to be used if `mask`
        is ``True``.
    verbose : bool, optional
        Whether to print message if colors are clipped.
    """
    # Notes:
    # I could use `numpy.clip` (`matplotlib.colors` uses this under the hood),
    # but we want to display messages. And anyway, premature efficiency is
    # the root of all evil, we're manipulating like 1000 colors max here, so
    # it's no big deal.
    message = 'Invalid' if mask else 'Clipped'
    colors = np.array(colors) # easier
    under = (colors<0)
    over  = (colors>1)
    if mask:
        colors[(under | over)] = gray
    else:
        colors[under] = 0
        colors[over]  = 1
    if verbose:
        for i,name in enumerate('rgb'):
            if under[:,i].any():
                warnings.warn(f'{message} "{name}" channel (<0).')
            if over[:,i].any():
                warnings.warn(f'{message} "{name}" channel (>1).')
    return colors

def _clip_cmap(cmap, left=None, right=None, N=None):
    """Helper function that cleanly divides linear segmented colormaps and
    subsamples listed colormaps. Full documentation is in `Colormap`. Note
    `N` is an alias for `right` for `ListedColormap` maps."""
    # Optionally clip edges or resample map.
    if left is None and right is None and (N is None or isinstance(cmap, mcolors.LinearSegmentedColormap)):
        return cmap
    # Listed colormaps, simply truncate colors in the list
    if isinstance(cmap, mcolors.ListedColormap):
        if N is not None:
            slicer = slice(None, N)
        else:
            slicer = slice(left, right)
        try:
            cmap = mcolors.ListedColormap(cmap.colors[slicer])
        except Exception:
            raise ValueError(f'Invalid indices {slicer} for listed colormap.')
        return cmap
    # Trickier procedure for segment data maps
    # First get segmentdata and parse input
    kwargs = {}
    data = cmap._segmentdata
    newdata = {}
    if left is None:
        left = 0
    if right is None:
        right = 1
    if hasattr(cmap, 'space'):
        kwargs['space'] = cmap.space
    # Next resample the segmentdata arrays
    dict_ = {key:value for key,value in data.items() if
            key not in ('gamma1','gamma2')}
    channels = {'saturation':'gamma1', 'luminance':'gamma2'}
    for key,xyy in dict_.items():
        # Get coordinates
        xyy     = np.array(xyy)
        x       = xyy[:,0]
        xleft,  = np.where(x>left)
        xright, = np.where(x<right)
        if len(xleft)==0:
            raise ValueError(f'Invalid x minimum {left}.')
        if len(xright)==0:
            raise ValueError(f'Invalid x maximum {right}.')
        # Slice
        # l is the first point where x>0 or x>left, should be at least 1
        # r is the last point where r<1 or r<right
        l, r = xleft[0], xright[-1]
        newxyy = xyy[l:r+1,:].copy()
        xl = xyy[l-1,1:] + (left - x[l-1])*(xyy[l,1:] - xyy[l-1,1:])/(x[l] - x[l-1])
        newxyy = np.concatenate(([[left, *xl]], newxyy), axis=0)
        xr = xyy[r,1:] + (right - x[r])*(xyy[r+1,1:] - xyy[r,1:])/(x[r+1] - x[r])
        newxyy = np.concatenate((newxyy, [[right, *xr]]), axis=0)
        newxyy[:,0] = (newxyy[:,0] - left)/(right - left)
        newdata[key] = newxyy
        # Retain the corresponding 'gamma' *segments*
        # Need more testing but so far so good
        if key in channels:
            name = channels[key]
            gamma = data[name]
            if np.iterable(gamma):
                gamma = gamma[l-1:r+1]
            newdata[name] = gamma
    # And finally rebuild map
    cmap = type(cmap)(cmap.name, newdata, **kwargs)
    return cmap

def _merge_cmaps(*args, ratios=1, name='merged', N=512, **kwargs):
    """Merges arbitrary colormaps. This is used when you pass multiple `args`
    to the `Colormap` function. Full documentation is in `Colormap`."""
    # Initial
    if len(args)<=1:
        raise ValueError('Need two or more input cmaps.')
    ratios = ratios or 1
    if isinstance(ratios, Number):
        ratios = [1]*len(args)
    imaps = [Colormap(cmap, N=None, **kwargs) for cmap in args] # set N=None to disable resamping

    # Simple process for listed colormap, just combine the colors
    if all(isinstance(cmap, mcolors.ListedColormap) for cmap in imaps):
        if not np.all(ratios==1):
            raise ValueError(f'Cannot assign different ratios when mering ListedColormaps.')
        colors = [color for cmap in imaps for color in cmap.colors]
        return mcolors.ListedColormap(colors, name=name, N=len(colors))
    elif not all(isinstance(cmap,mcolors.LinearSegmentedColormap) for cmap in imaps):
        raise ValueError('All colormaps must be of the same type (Listed or LinearSegmented).')

    # More complex for continuous maps
    # Get the segmentdata, and figure out x-coordinates for merging
    keys = {key for cmap in imaps for key in cmap._segmentdata.keys()
            if key not in ('gamma1','gamma2')}
    ratios = np.array(ratios)/np.sum(ratios) # so if 4 cmaps, will be 1/4
    x0 = np.concatenate([[0], np.cumsum(ratios)])
    xw = x0[1:] - x0[:-1] # weights for averages
    # Mixed types
    types = {type(cmap) for cmap in imaps}
    if len(types)>1:
        raise ValueError(f'Mixed colormap types {types}. Maps must all be LinearSegmentedColormap or PerceptuallyUniformColormap.')
    type_ = types.pop()

    # Special consideration for PerceptuallyUniformColormaps
    kwargs = {}
    segmentdata = {}
    if type_ is PerceptuallyUniformColormap:
        # Combine gamma values, where gamma1 corresponds to saturation and
        # gamma2 to luminance. Account for scalar gamma or gamma list options.
        channels = {0:'saturation', 1:'luminance'}
        for i,key in enumerate(('gamma1', 'gamma2')):
            if key not in segmentdata:
                segmentdata[key] = []
            for cmap in imaps:
                gamma = cmap._segmentdata[key]
                if not np.iterable(gamma):
                    gamma = [gamma]*(len(cmap._segmentdata[channels[i]])-1) # length is *number* of rows in segmentdata
                segmentdata[key].extend([*gamma])
        # Get the colorspace
        spaces = {cmap.space for cmap in imaps}
        if len(spaces)>1:
            raise ValueError(f'Cannot merge colormaps in the different HSL spaces {repr(spaces)}.')
        kwargs['space'] = spaces.pop()
    # Combine the segmentdata, and use the y1/y2 slots at merge points
    # so the transition is immediate (i.e. the colormap will never try to
    # interpolate between the end colors of different colormaps)
    for key in keys:
        # Combine xyy data
        datas = []
        callable_ = [callable(cmap._segmentdata[key]) for cmap in imaps]
        if not all(callable_) and any(callable_):
            raise ValueError('Mixed callable and non-callable colormap values.')
        if all(callable_): # expand range from x-to-w to 0-1
            for x,w,cmap in zip(x0[:-1], xw, imaps):
                data = lambda x: data((x - x0)/w) # WARNING: untested!
                datas.append(data)
        else:
            for x,w,cmap in zip(x0[:-1], xw, imaps):
                data = np.array(cmap._segmentdata[key])
                data[:,0] = x + w*data[:,0]
                datas.append(data)
            for i in range(len(datas)-1):
                datas[i][-1,2] = datas[i+1][0,2]
                datas[i+1] = datas[i+1][1:,:]
            data = np.concatenate(datas, axis=0)
            data[:,0] = data[:,0]/data[:,0].max(axis=0) # scale to make maximum exactly 1 (avoid floating point errors)
        segmentdata[key] = data
    return type(imaps[0])(name, segmentdata, N=N, **kwargs)

def _make_segmentdata_array(values, ratios=None, reverse=False, **kwargs):
    """Constructs a list of linear segments for an individual channel.
    This was made so that user can input e.g. a callable function for
    one channel, but request linear interpolation for another one."""
    # Handle function handles
    if callable(values):
        if reverse:
            values = lambda x: values(1-x)
        return values # just return the callable
    values = np.atleast_1d(values)
    if len(values)==1:
        value = values[0]
        return [(0, value, value), (1, value, value)] # just return a constant transition

    # Get x coordinates
    if not np.iterable(values):
        raise TypeError('Colors must be iterable.')
    if ratios is not None:
        xvals = np.atleast_1d(ratios) # could be ratios=1, i.e. dummy
        if len(xvals) != len(values) - 1:
            raise ValueError(f'Got {len(values)} values, but {len(ratios)} ratios.')
        xvals = np.concatenate(([0], np.cumsum(xvals)))
        xvals = xvals/np.max(xvals) # normalize to 0-1
    else:
        xvals = np.linspace(0,1,len(values))

    # Build vector
    array = []
    slicer = slice(None,None,-1) if reverse else slice(None)
    for x,value in zip(xvals,values[slicer]):
        array.append((x, value, value))
    return array

def make_mapping_array(N, data, gamma=1.0, reverse=False):
    r"""
    Mostly a copy of `~matplotlib.colors.makeMappingArray`, but allows
    *circular* hue gradations along 0-360, disables clipping of
    out-of-bounds channel values, and with fancier "gamma" scaling.

    Parameters
    ----------
    N : int
        Number of points in the generated lookup table.
    data : 2D array-like
        List of :math:`(x, y_0, y_1)` tuples specifying the channel jump (from
        :math:`y_0` to :math:`y_1`) and the :math:`x` coordinate of that
        transition (ranges between 0 and 1).
        See `~matplotlib.colors.LinearSegmentedColormap` for details.
    gamma : float or list of float, optional
        To obtain channel values between coordinates :math:`x_i` and
        :math:`x_{i+1}` in rows :math:`i` and :math:`i+1` of `data`,
        we use the formula:

        .. math::

            y = y_{1,i} + w_i^{\gamma_i}*(y_{0,i+1} - y_{1,i})

        where :math:`\gamma_i` corresponds to `gamma` and the weight
        :math:`w_i` ranges from 0 to 1 between rows ``i`` and ``i+1``.
        If `gamma` is float, it applies to every transition. Otherwise,
        its length must equal ``data.shape[0]-1``.
    reverse : bool, optional
        If ``True``, :math:`w_i^{\gamma_i}` is replaced with
        :math:`1 - (1 - w_i)^{\gamma_i}` -- that is, when `gamma` is greater
        than 1, this weights colors toward *higher* channel values instead
        of lower channel values.

        This is implemented in case we want to apply *equal* "gamma scaling"
        to different HSL channels in different directions. Usually, this
        is done to weight low data values with higher luminance *and* lower
        saturation, thereby emphasizing "extreme" data values with stronger
        colors.
    """
    # Optionally allow for ***callable*** instead of linearly interpolating
    # between line segments
    gammas = np.atleast_1d(gamma)
    if (gammas < 0.01).any() or (gammas > 10).any():
        raise ValueError('Gamma can only be in range [0.01,10].')
    if callable(data):
        if len(gammas)>1:
            raise ValueError('Only one gamma allowed for functional segmentdata.')
        x = np.linspace(0, 1, N)**gamma
        lut = np.array(data(x), dtype=float)
        return lut

    # Get array
    try:
        data = np.array(data)
    except Exception:
        raise TypeError('Data must be convertible to an array.')
    shape = data.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError('Data must be nx3 format.')
    if len(gammas)!=1 and len(gammas)!=shape[0]-1:
        raise ValueError(f'Need {shape[0]-1} gammas for {shape[0]}-level mapping array, but got {len(gamma)}.')
    if len(gammas)==1:
        gammas = np.repeat(gammas, shape[:1])

    # Get indices
    x  = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    if x[0] != 0.0 or x[-1] != 1.0:
        raise ValueError('Data mapping points must start with x=0 and end with x=1.')
    if (np.diff(x) < 0).any():
        raise ValueError('Data mapping points must have x in increasing order.')
    x = x*(N - 1)

    # Get distances from the segmentdata entry to the *left* for each requested
    # level, excluding ends at (0,1), which must exactly match segmentdata ends
    xq = (N - 1)*np.linspace(0, 1, N)
    ind = np.searchsorted(x, xq)[1:-1] # where xq[i] must be inserted so it is larger than x[ind[i]-1] but smaller than x[ind[i]]
    distance = (xq[1:-1] - x[ind - 1])/(x[ind] - x[ind - 1])
    # Scale distances in each segment by input gamma
    # The ui are starting-points, the ci are counts from that point
    # over which segment applies (i.e. where to apply the gamma)
    _, uind, cind = np.unique(ind, return_index=True, return_counts=True)
    for i,(ui,ci) in enumerate(zip(uind,cind)): # i will range from 0 to N-2
        # Test if 1
        gamma = gammas[ind[ui]-1] # the relevant segment is to *left* of this number
        if gamma==1:
            continue
        # By default, weight toward a *lower* channel value (i.e. bigger
        # exponent implies more colors at lower value)
        # Again, the relevant 'segment' is to the *left* of index returned by searchsorted
        ir = False
        if ci>1: # i.e. more than 1 color in this 'segment'
            ir = ((y0[ind[ui]] - y1[ind[ui]-1]) < 0) # by default want to weight toward a *lower* channel value
        if reverse:
            ir = (not ir)
        if ir:
            distance[ui:ui + ci] = 1 - (1 - distance[ui:ui + ci])**gamma
        else:
            distance[ui:ui + ci] **= gamma

    # Perform successive linear interpolations all rolled up into one equation
    lut = np.zeros((N,), float)
    lut[1:-1] = distance*(y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0]  = y1[0]
    lut[-1] = y0[-1]
    return lut

#------------------------------------------------------------------------------#
# Generalized colormap/cycle constructors
#------------------------------------------------------------------------------#
def colors(*args, **kwargs):
    """Alias for `Cycle`."""
    return Cycle(*args, **kwargs)

def Colormap(*args, name=None, cyclic=False, N=None,
        cut=None, left=None, right=None, x=None, reverse=False,
        ratios=1, gamma=None, gamma1=None, gamma2=None,
        save=False,
        **kwargs):
    """
    Convenience function for generating and **merging** colormaps
    in a variety of ways.

    Parameters
    ----------
    *args : `~matplotlib.colors.Colormap`, dict, list of str, or str
        Each arg generates a single colormap. If ``len(args)>1``, the colormaps
        are merged.

        If the arg is a `~matplotlib.colors.Colormap`, nothing more is done.
        Otherwise, the colormap is generated as follows:

        * If arg is a str and is a "registered" colormap name, that colormap
          is used.
        * If arg is a str and is a "registered" color name, a
          monochromatic colormap is generated with `monochrome_cmap`.
        * If arg is a list of str, it is assumed a list of color names or
          hex names, and is used to make a `~matplotlib.colors.ListedColormap`.
        * If arg is dict, there are two options: if the dict contains the keys
          ``'red'``, ``'green'``, and ``'blue'``, it is passed to the
          `~matplotlib.colors.LinearSegmentedColormap` initializer. Otherwise,
          the dict is passed to the `~PerceptuallyUniformColormap.from_hsl`
          `PerceptuallyUniformColormap` constructor.

        For the monochromatic colormaps, the color name string can also
        look like ``'name90'``, where the trailing
        number indicates the maximum luminance of the colormap. By default,
        this is 100 (i.e. pure white).
    name : None or str, optional
        Name of colormap. Default name is ``'no_name'``.
        The resulting colormap can then be invoked by passing ``cmap='name'``
        to plotting functions like `~matplotlib.axes.Axes.contourf`.
    cyclic : bool, optional
        Whether the colormap is cyclic. Will cause `~proplot.axes.wrapper_cmap`
        to pass this flag to `BinNorm`. This will prevent having the same color
        on either end of the colormap.
    N : None or int, optional
        Number of colors to generate in the hidden lookupt table ``_lut``.
        By default, a relatively high resolution of 256 is chosen (see notes).
    cut : None or float, optional
        Optionally cut out colors in the **center** of the colormap. This is
        extremely useful for diverging colormaps, in case you want to have a
        sharper cutoff between negative and positive values. For example,
        ``cut=0.1`` cuts out the middle 10% of the colormap.
    left, right : None or float or list of float, optional
        Optionally *delete* colors on the left and right sides of the
        colormap(s). For example, ``left=0.1`` deletes the leftmost
        10% of the colormap; ``right=0.9`` deletes the rightmost 10%.

        If list, length must match ``len(args)``, and applies to *each*
        colormap in the list before they are merged. If float, applies
        to the *final* colormap. No difference if ``len(args)`` is 1.
    x : None or (float, float), optional
        List containing the keyword args `(left, right)`.
    reverse : bool or list of bool, optional
        Optionally reverse the colormap(s).

        If list, length must match ``len(args)``, and applies to *each*
        colormap in the list before they are merged. If bool, applies
        to the *final* colormap. No difference if ``len(args)`` is 1.
    ratios : list of float, optional
        Indicates the ratios used to *combine* the colormaps. Length must
        equal ``len(args)`` (ignored if ``len(args)`` is 1).

        For example, if `args` contains ``['blues', 'reds']`` and `ratios`
        is ``[2, 1]``, this generates a colormap with two-thirds blue
        colors on the left and one-third red colors in the right.
    gamma1, gamma2, gamma : float, optional
        Gamma-scaling for the saturation, luminance, and both channels
        for perceptualy uniform colormaps. See the
        `PerceptuallyUniformColormap` documentation.

        So far there is no way to apply this differently for each colormap
        in `args`.
    save : bool, optional
        Whether to save the colormap in the folder ``~/.proplot``. The
        folder is created if it does not already exist.

        If the colormap is a `~matplotlib.colors.ListedColormap` (i.e. a
        "color cycle"), the list of hex strings are written to ``name.hex``.

        If the colormap is a `~matplotlib.colors.LinearSegmentedColormap`,
        the segment data dictionary is written to ``name.json``.

    Note
    ----
    The resampling method described in `this post
    <https://stackoverflow.com/q/48613920/4970632>`_ is not great. All it does is
    reduce the lookup table size -- what ends up happening under the hood is
    matplotlib tries to *evenly* draw ``N-1`` (``'min'`` or ``'max'``) or ``N-2``
    (``'neither'``) colors from a lookup table with ``N`` colors, which means
    it simply *skips over* 1 or 2 colors in the middle of the lookup table,
    which will cause visual jumps!

    Instead, we make the stored lookup table completely independent from the
    number of levels. It will have many high-res segments while the colormap
    `N` is very small.
    """
    # Initial stuff
    N_ = N or rcParams['image.lut']
    imaps = []
    name = name or 'no_name' # must have name, mcolors utilities expect this
    if x is not None:
        left, right = x
    if len(args)==0:
        args = [rcParams['image.cmap']] # use default
    for i,cmap in enumerate(args):
        # Retrieve Colormap instance. Makes sure lookup table is reset.
        if cmap is None:
            cmap = rcParams['image.cmap']
        if isinstance(cmap,str) and cmap in mcm.cmap_d:
            cmap = mcm.cmap_d[cmap]
        if isinstance(cmap, mcolors.LinearSegmentedColormap):
            # Resample, allow overriding the gamma
            # Copy over this add-on attribute
            # NOTE: Calling resample means cmaps are un-initialized
            cyclic = getattr(cmap, '_cyclic', False)
            cmap = cmap._resample(N_)
            cmap._cyclic = cyclic
            if isinstance(cmap, PerceptuallyUniformColormap):
                gamma1 = _default(gamma, gamma1)
                gamma2 = _default(gamma, gamma2)
                segmentdata = cmap._segmentdata
                if gamma1:
                    segmentdata['gamma1'] = gamma1
                if gamma2:
                    segmentdata['gamma2'] = gamma2
            elif gamma:
                cmap._gamma = gamma
        elif isinstance(cmap, mcolors.ListedColormap):
            pass
        # Build colormap on-the-fly
        elif isinstance(cmap, dict):
            # Dictionary of hue/sat/luminance values or 2-tuples representing linear transition
            if {*cmap.keys()} == {'red','green','blue'}:
                cmap = mcolors.LinearSegmentedColormap(name, cmap, N=N_)
            else:
                cmap = PerceptuallyUniformColormap.from_hsl(name, N=N_, **cmap)
        elif not isinstance(cmap, str) and np.iterable(cmap) and all(np.iterable(color) for color in cmap):
            # List of color tuples or color strings, i.e. iterable of iterables
            # Transform C0, C1, etc. to their actual names first
            cmap = [_transform_cycle(color) for color in cmap]
            cmap = mcolors.ListedColormap(cmap, name=name, **kwargs)
        else:
            # Monochrome colormap based from input color (i.e. single hue)
            # TODO: What if colormap names conflict with color names! Maybe address
            # this! Currently, this makes it impossible to make a monochrome colormap
            # from some named color if that name also exists for a colormap.
            # Try to convert to RGB
            cmap = _transform_cycle(cmap)
            fade = kwargs.pop('fade', 90) 
            if isinstance(cmap, str): # not a color tuple
                regex = '([0-9]+)$'
                match = re.search(regex, cmap) # declare maximum luminance with e.g. red90, blue70, etc.
                cmap = re.sub(regex, '', cmap) # remove options
                if match:
                    fade = float(match.group(1)) # default fade to 100 luminance
            # Build colormap
            cmap = to_rgb(cmap) # to ensure is hex code/registered color
            cmap = monochrome_cmap(cmap, fade, name=name, N=N_, **kwargs)

        # Optionally transform colormap by clipping colors or reversing
        if np.iterable(reverse) and reverse[i]:
            cmap = cmap.reversed()
        cmap = _clip_cmap(cmap, None if not np.iterable(left) else left[i],
                                None if not np.iterable(right) else right[i], N=N)
        imaps += [cmap]

    # Now merge the result of this arbitrary user input
    # Since we are merging cmaps, potentially *many* color transitions; use big number by default
    N_ = N_*len(imaps)
    if len(imaps)>1:
        cmap = _merge_cmaps(*imaps, name=name, ratios=ratios, N=N_)

    # Transform colormap
    if not np.iterable(reverse) and reverse:
        cmap = cmap.reversed()
    left = None if np.iterable(left) else left
    right = None if np.iterable(right) else right
    # Cut out middle colors of a diverging map
    if cut: # non-zero and not None
        cright, cleft = 0.5 - cut/2, 0.5 + cut/2
        lcmap = _clip_cmap(cmap, left, cright)
        rcmap = _clip_cmap(cmap, cleft, right)
        cmap = _merge_cmaps(lcmap, rcmap, name=name, N=N_)
    # Cut out either edge
    else:
        cmap = _clip_cmap(cmap, left, right, N=N)
    # Add property
    # Only a couple builtin maps should already have this attribute, from register_cmaps
    if not hasattr(cmap, '_cyclic'):
        cmap._cyclic = cyclic
    # Initialize (the _resample methods generate new colormaps,
    # so current one is uninitializied)
    if not cmap._isinit:
        cmap._init()

    # Perform crude resampling of data, i.e. just generate a low-resolution
    # lookup table instead.
    # NOTE: This approach is no longer favored; instead we generate hi-res
    # lookup table and use BinNorm to discretize colors, much more flexible.
    # if isinstance(cmap, mcolors.LinearSegmentedColormap) and N is not None:
    #     offset = {'neither':-1, 'max':0, 'min':0, 'both':1}
    #     if extend not in offset:
    #         raise ValueError(f'Unknown extend option {extend}.')
    #     cmap = cmap._resample(N - offset[extend]) # see mcm.get_cmap source

    # Register the colormap
    mcm.cmap_d[name] = cmap
    # Optionally save colormap to disk
    if save:
        if not os.path.isdir(_data_user):
            os.mkdir(_data_user)
        # Save listed colormap i.e. color cycle
        if isinstance(cmap, mcolors.ListedColormap):
            basename = f'{name}.hex'
            filename = os.path.join(_data_user, basename)
            with open(filename, 'w') as f:
                f.write(','.join(mcolors.to_hex(color) for color in cmap.colors))
        # Save segment data directly
        else:
            basename = f'{name}.json'
            filename = os.path.join(_data_user, basename)
            data = {}
            for key,value in cmap._segmentdata.items():
                data[key] = value.astype(float).tolist() # from np.float to builtin float, and to list of lists
            if hasattr(cmap, 'space'):
                data['space'] = cmap.space
            with open(filename, 'w') as file:
                json.dump(data, file, indent=4)
        print(f'Saved colormap to "{basename}".')
    return cmap

def Cycle(*args, samples=10, vmin=0, vmax=1, **kwargs):
    """
    Convenience function that builds lists of colors from colormaps or returns
    lists of colors from existing registered cycles.

    For color cycles (i.e. `~matplotlib.colors.ListedColormap` instances),
    we just select colors from that list.
    For colormaps (i.e. `~matplotlib.colors.LinearSegmentedColormap` instances),
    we draw samples from the full range of colormap colors.

    Parameters
    ----------
    *args
        Passed to `Colormap`. If ``args[-1]`` is a float or list of float,
        it is used as the `samples` argument.

        This allows the user to declare new color cycles with, for example,
        ``proplot.rc.cycle = ['blues', 'reds', 20]``
    samples : float or list of float, optional
        Last element of `args`; detected automatically if this element is
        iterable or a number.

        This implementation seems weird, but this allows the user to do
        stuff like ``proplot.rc.cycle = ['blues', 'reds', 20]``
        Required if `Colormap` returns a smooth
        colormap (i.e. a `~matplotlib.colors.LinearSegmentedColormap` 
        instance).
    vmin, vmax : float, optional
        The minimum and maximum data values, used to scale `samples`.

    Other parameters
    ----------------
    **kwargs
        Passed to `Colormap`.

    In the latter case, we will draw samples from that colormap by (default)
    drawing from Use vmin/vmax to scale your samples.
    """
    # Two modes:
    # 1) User inputs some number of samples; 99% of time, use this
    # to get samples from a LinearSegmentedColormap draw colors.
    if isinstance(args[-1], Number) or \
            (np.iterable(args[-1]) and not isinstance(args[-1], (str, dict))):
        args, samples = args[:-1], args[-1]
    # 2) User inputs a simple list; 99% of time, use this
    # to build up a simple ListedColormap.
    elif len(args)>1:
        args = [args] # presumably send a list of colors
    cmap = Colormap(*args, **kwargs) # the cmap object itself
    name = cmap.name
    if isinstance(cmap, mcolors.ListedColormap):
        # Just get the colors
        colors = cmap.colors
    elif isinstance(cmap, mcolors.LinearSegmentedColormap): # or subclass
        # Employ ***more flexible*** version of get_cmap() method, which does this:
        # LinearSegmentedColormap(self.name, self._segmentdata, lutsize)
        if isinstance(samples, Number):
            samples = np.linspace(0, 1, samples) # from edge to edge
        elif np.iterable(samples):
            samples = np.array(samples)
        else:
            raise ValueError(f'Invalid samples "{samples}".')
        colors = cmap((samples-vmin)/(vmax-vmin))
    else:
        raise ValueError(f'Colormap returned weird object type: {type(cmap)}.')
    # Return
    colors = [tuple(color) if not isinstance(color,str) else color for color in colors]
    return CycleList(colors, name)

class PerceptuallyUniformColormap(mcolors.LinearSegmentedColormap):
    """
    Similar to `~matplotlib.colors.LinearSegmentedColormap`, but instead
    of varying the RGB channels, we vary hue, saturation, and luminance in
    either the perceptually uniform HCL colorspace or the HSLuv or HPLuv
    scalings of HCL.
    """
    def __init__(self, name, segmentdata, space='hsl', mask=False,
        gamma=None, gamma1=None, gamma2=None, **kwargs):
        """
        Parameters
        ----------
        name : str
            The colormap name.
        segmentdata : dict-like
            Dictionary mapping containing the keys ``'hue'``, ``'saturation'``,
            and ``'luminance'``. Values should be lists containing any of
            the following channel specifiers:

                1. Numbers, within the range 0-360 for hue and 0-100 for
                   saturation and luminance.
                2. Color string names or hex tags, in which case the channel
                   value for that color is looked up.

            See `~matplotlib.colors.LinearSegmentedColormap` for details.
        space : {'hsl', 'hcl', 'hpl'}, optional
            The hue, saturation, luminance-style colorspace to use for
            interpreting the channels. See `this page
            <http://www.hsluv.org/comparison/>`_ for a description of each
            colorspace.
        mask : bool, optional
            When we interpolate across HSL space, we can end
            up with "impossible" RGB colors (colors with channel values
            >1).

            If `mask` is ``True``, these "impossible" colors are masked
            out as black. Otherwise, the channels are just clipped to 1.
            Default is ``False``.
        gamma1 : None or float, optional
            If >1, makes low saturation colors more prominent. If <1,
            makes high saturation colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.
        gamma2 : None or float, optional
            If >1, makes high luminance colors more prominent. If <1,
            makes low luminance colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.
        gamma : None or float, optional
            Use this to identically set `gamma1` and `gamma2` at once.

        Example
        -------
        The following is a valid `segmentdata` dictionary, using color string
        names for the hue instead of numbers between 0 and 360.

        .. code-block:: python

            dict(hue       = [[0, 'red', 'red'], [1, 'blue', 'blue']],
                saturation = [[0, 100, 100], [1, 100, 100]],
                luminance  = [[0, 100, 100], [1, 20, 20]])

        Note
        ----
        `gamma1` emphasizes *low* saturation colors and `gamma2` emphasizes
        *high* luminance colors because this seems to be what the
        ColorBrewer2.0 maps do.  "White" and "pale" values at the center of
        diverging maps and on the left of sequential maps are given
        extra emphasis.
        """
        # Attributes
        # NOTE: Don't allow power scaling for hue because that would be weird.
        # Idea is want to allow skewing so dark/saturated colors are
        # more isolated/have greater intensity.
        # NOTE: We add gammas to the segmentdata dictionary so it can be
        # pickled into .npy file
        space = _get_space(space)
        if 'gamma' in kwargs:
            raise ValueError('Standard gamma scaling disabled. Use gamma1 or gamma2 instead.')
        gamma1 = _default(gamma, gamma1)
        gamma2 = _default(gamma, gamma2)
        segmentdata['gamma1'] = _default(gamma1, segmentdata.get('gamma1', None), 1.0)
        segmentdata['gamma2'] = _default(gamma2, segmentdata.get('gamma2', None), 1.0)
        self._space = space
        self._mask  = mask
        # First sanitize the segmentdata by converting color strings to their
        # corresponding channel values
        keys   = {*segmentdata.keys()}
        target = {'hue', 'saturation', 'luminance', 'gamma1', 'gamma2'}
        if keys != target and keys != {*target, 'alpha'}:
            raise ValueError(f'Invalid segmentdata dictionary with keys {keys}.')
        for key,array in segmentdata.items():
            # Allow specification of channels using registered string color names
            if 'gamma' in key:
                continue
            if callable(array):
                continue
            for i,xyy in enumerate(array):
                xyy = list(xyy) # make copy!
                for j,y in enumerate(xyy[1:]): # modify the y values
                    xyy[j+1] = _get_channel(y, key, space)
                segmentdata[key][i] = xyy
        # Initialize
        # NOTE: Our gamma1 and gamma2 scaling is just fancy per-channel
        # gamma scaling, so disable the standard version.
        super().__init__(name, segmentdata, gamma=1.0, **kwargs)

    def reversed(self, name=None):
        """Returns reversed colormap."""
        if name is None:
            name = self.name + '_r'
        def factory(dat):
            def func_r(x):
                return dat(1.0 - x)
            return func_r
        data_r = {}
        for key,xyy in self._segmentdata.items():
            if key in ('gamma1', 'gamma2', 'space'):
                if 'gamma' in key: # optional per-segment gamma
                    xyy = np.atleast_1d(xyy)[::-1]
                data_r[key] = xyy
                continue
            elif callable(xyy):
                data_r[key] = factory(xyy)
            else:
                data_r[key] = [[1.0 - x, y1, y0] for x, y0, y1 in reversed(xyy)]
        return PerceptuallyUniformColormap(name, data_r, space=self._space)

    def _init(self):
        """As with `~matplotlib.colors.LinearSegmentedColormap`, but converts
        each value in the lookup table from 'input' to RGB."""
        # First generate the lookup table
        channels = ('hue','saturation','luminance')
        reverse = (False, False, True) # gamma weights *low chroma* and *high luminance*
        gammas = (1.0, self._segmentdata['gamma1'], self._segmentdata['gamma2'])
        self._lut_hsl = np.ones((self.N+3, 4), float) # fill
        for i,(channel,gamma,reverse) in enumerate(zip(channels, gammas, reverse)):
            self._lut_hsl[:-3,i] = make_mapping_array(self.N, self._segmentdata[channel], gamma, reverse)
        if 'alpha' in self._segmentdata:
            self._lut_hsl[:-3,3] = make_mapping_array(self.N, self._segmentdata['alpha'])
        self._lut_hsl[:-3,0] %= 360
        # self._lut_hsl[:-3,0] %= 359 # wrong
        # Make hues circular, set extremes (i.e. copy HSL values)
        self._lut = self._lut_hsl.copy() # preserve this, might want to check it out
        self._set_extremes() # generally just used end values in segmentdata
        self._isinit = True
        # Now convert values to RGBA, and clip colors
        for i in range(self.N+3):
            self._lut[i,:3] = to_rgb(self._lut[i,:3], self._space)
        self._lut[:,:3] = _clip_colors(self._lut[:,:3], self._mask)

    def _resample(self, N):
        """Returns a new colormap with *N* entries."""
        return PerceptuallyUniformColormap(self.name, self._segmentdata, self._space, self._mask, N=N)

    @staticmethod
    def from_hsl(name, h=0, s=100, l=[100, 20], c=None, a=None,
            hue=None, saturation=None, luminance=None, chroma=None, alpha=None,
            ratios=None, reverse=False, **kwargs):
        """
        Makes a `~PerceptuallyUniformColormap` by specifying the hue, saturation,
        and luminance transitions individually.

        Parameters
        ----------
        h, hue : float, str, or list thereof, optional
            Hue channel value or list of values. Values can be
            any of the following:

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex tags, in which case the channel
               value for that color is looked up.

            If scalar, the hue does not change across the colormap.
        s, saturation : float, str, or list thereof, optional
            As with hue, but for the saturation channel.
        l, luminance, c, chroma : float, str, or list thereof, optional
            As with hue, but for the luminance channel.
        a, alpha : float, str, or list thereof, optional
            As with hue, but for the alpha channel (the transparency).
        ratios : None or list of float, optional
            Relative extent of the transitions indicated by the channel
            value lists.

            For example, ``luminance=[100,50,0]`` with ``ratios=[2,1]``
            places the *x*-coordinate where the luminance is 50 at 0.66 --
            the white to gray transition is "slower" than the gray to black
            transition.
        reverse : bool, optional
            Whether to reverse the final colormap.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.

        Todo
        ----
        Add ability to specify *discrete "jump" transitions*. Currently, this
        only lets you generate colormaps with smooth transitions, and is not
        necessarily suited for e.g. for making diverging colormaps.
        """
        # Build dictionary, easy peasy
        h = _default(hue, h)
        s = _default(chroma, c, saturation, s)
        l = _default(luminance, l)
        a = _default(alpha, a, 1.0)
        cs = ['hue', 'saturation', 'luminance', 'alpha']
        channels = [h, s, l, a]
        cdict = {}
        for c,channel in zip(cs,channels):
            cdict[c] = _make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

    @staticmethod
    def from_list(name, color_list,
        ratios=None, reverse=False,
        **kwargs):
        """
        Makes a `PerceptuallyUniformColormap` from a list of (hue, saturation,
        luminance) tuples.

        Parameters
        ----------
        name : str
            The colormap name.
        color_list : list of length-3 tuples
            List containing HSL color tuples. The tuples can contain any
            of the following channel value specifiers:

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex tags, in which case the channel
               value for that color is looked up.

        ratios : None or list of float, optional
            Length ``len(color_list)-1`` list of scales for *x*-coordinate
            transitions between colors. Bigger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.
        reverse : bool, optional
            Whether to reverse the result.
        """
        # Dictionary
        cdict = {}
        channels = [*zip(*color_list)]
        if len(channels) not in (3,4):
            raise ValueError(f'Bad color list: {color_list}')
        cs = ['hue', 'saturation', 'luminance']
        if len(channels)==4:
            cs += ['alpha']
        else:
            cdict['alpha'] = 1.0 # dummy function that always returns 1.0
        # Build data arrays
        for c,channel in zip(cs,channels):
            cdict[c] = _make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

def monochrome_cmap(color, fade, reverse=False, space='hpl', name='monochrome', **kwargs):
    """
    Makes a monochromatic "sequential" colormap that blends from near-white
    to the input color.

    Parameters
    ----------
    color : str or (R,G,B) tuple
        Color RGB tuple, hex string, or named color string.
    fade : float or str or (R,G,B) tuple
        The luminance channel strength, or color from which to take the luminance channel.
    reverse : bool, optional
        Whether to reverse the colormap.
    space : {'hsl', 'hcl', 'hpl'}, optional
        Colorspace in which the luminance is varied.
    name : str, optional
        Colormap name. Default is ``'monochrome'``.

    Other parameters
    ----------------
    **kwargs
        Passed to `PerceptuallyUniformColormap.from_hsl` static method.
    """
    # Get colorspace
    # NOTE: If you use HSL space, will get very saturated colors in the middle
    # of the map (around 50% luminance); otherwise chroma won't change
    h, s, l = to_xyz(to_rgb(color), space)
    if isinstance(fade, Number): # allow just specifying the luminance channel
        fs, fl = s, fade # fade to *same* saturation by default
        fs = s/2
    else:
        _, fs, fl = to_xyz(to_rgb(fade), space)
    index = slice(None,None,-1) if reverse else slice(None)
    return PerceptuallyUniformColormap.from_hsl(name, h,
            [fs,s][index], [fl,l][index], space=space, **kwargs)

#------------------------------------------------------------------------------#
# Return arbitrary normalizer
#------------------------------------------------------------------------------
def Norm(norm_in, levels=None, values=None, norm=None, **kwargs):
    """
    Returns an arbitrary `~matplotlib.colors.Normalize` instance.

    Parameters
    ----------
    norm_in : str or `~matplotlib.colors.Normalize`
        Key name for the normalizer. The recognized normalizer key names
        are as follows:

        ===============  =================================
        Key              Class
        ===============  =================================
        ``'none'``       `~matplotlib.colors.NoNorm`
        ``'null'``       `~matplotlib.colors.NoNorm`
        ``'zero'``       `MidpointNorm`
        ``'midpoint'``   `MidpointNorm`
        ``'segments'``   `LinearSegmentedNorm`
        ``'segmented'``  `LinearSegmentedNorm`
        ``'log'``        `~matplotlib.colors.LogNorm`
        ``'linear'``     `~matplotlib.colors.Normalize`
        ``'power'``      `~matplotlib.colors.PowerNorm`
        ``'symlog'``     `~matplotlib.colors.SymLogNorm`
        ===============  =================================

    levels, values : array-like
        The level edges (`levels`) or centers (`values`) passed
        to `LinearSegmentedNorm`.
    norm : None or normalizer spec, optional
        The normalizer for *pre-processing*. Used only for
        the `LinearSegmentedNorm` normalizer.

    Other parameters
    ----------------
    **kwargs
        Passed to the `~matplotlib.colors.Normalize` initializer.
        See `this tutorial <https://matplotlib.org/tutorials/colors/colormapnorms.html>`_
        for more info.
    """
    norm, norm_preprocess = norm_in, norm
    if isinstance(norm, mcolors.Normalize):
        return norm
    if levels is None and values is not None:
        levels = utils.edges(values)
    if isinstance(norm, str):
        # Get class
        norm_out = normalizers.get(norm, None)
        if norm_out is None:
            raise ValueError(f'Unknown normalizer "{norm}". Options are {", ".join(normalizers.keys())}.')
        # Instantiate class
        if norm_out is MidpointNorm:
            if not np.iterable(levels):
                raise ValueError(f'Need levels for normalizer "{norm}". Received levels={levels}.')
            kwargs.update({'vmin':min(levels), 'vmax':max(levels)})
        elif norm_out is LinearSegmentedNorm:
            if not np.iterable(levels):
                raise ValueError(f'Need levels for normalizer "{norm}". Received levels={levels}.')
            kwargs.update({'levels':levels, 'norm':norm_preprocess})
        norm_out = norm_out(**kwargs) # initialize
    else:
        raise ValueError(f'Unknown norm "{norm_out}".')
    return norm_out

#------------------------------------------------------------------------------
# Very important normalization class. Essentially there are two ways to create
# discretized color levels from a functional/segmented colormap:
#   1) Make lo-res lookup table.
#   2) Make hi-res lookup table, but discretize the lookup table indices
#      generated by your normalizer. This is what BinNorm does.
# Have found the second method was easier to implement/more flexible. So the
# below is used to always discretize colors.
#------------------------------------------------------------------------------
# WARNING: Many methods in ColorBarBase tests for class membership, crucially
# including _process_values(), which if it doesn't detect BoundaryNorm will
# end up trying to infer boundaries from inverse() method. So make it parent class.
class BinNorm(mcolors.BoundaryNorm):
    """
    This normalizer is used for all colormap plots. It can be thought of as a
    "parent" normalizer: it first scales the data according to any
    arbitrary `~matplotlib.colors.Normalize` class, then maps the normalized
    values ranging from 0-1 into **discrete** levels.

    This maps to colors by the closest **index** in the color list. Even if
    your levels edges are weirdly spaced (e.g. [-1000, 100, 0,
    100, 1000] or [0, 10, 12, 20, 22]), the "colormap coordinates" for these
    levels will be [0, 0.25, 0.5, 0.75, 1].

    Note
    ----
    If you are using a diverging colormap with ``extend='max'`` or
    ``extend='min'``, the center will get messed up. But that is very strange
    usage anyway... so please just don't do that :)
    """
    def __init__(self, levels, norm=None, clip=False, step=1.0, extend='neither'):
        """
        Parameters
        ----------
        levels : list of float
            The discrete data levels.
        norm : None or `~matplotlib.colors.Normalize`, optional
            The normalizer used to transform `levels` and all data passed
            to `BinNorm.__call__` *before* discretization.
        step : float, optional
            The intensity of the transition to out-of-bounds color, as a
            faction of the *average* step between in-bounds colors. The
            default is ``1``.
        extend : {'neither', 'both', 'min', 'max'}, optional
            Which direction colors will be extended. No matter the `extend`
            option, `BinNorm` ensures colors always extend through the
            extreme end colors.
        clip : bool, optional
            A `~matplotlib.colors.Normalize` option.
        """
        # Declare boundaries, vmin, vmax in True coordinates.
        # Notes:
        # * Idea is that we bin data into len(levels) discrete x-coordinates,
        #   and optionally make out-of-bounds colors the same or different
        # * Don't need to call parent __init__, this is own implementation
        #   Do need it to subclass BoundaryNorm, so ColorbarBase will detect it
        #   See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        levels = np.atleast_1d(levels)
        if levels.size<=1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError(f'Levels {levels} passed to Normalize() must be monotonically increasing.')
        if extend not in ('both','min','max','neither'):
            raise ValueError(f'Unknown extend option "{extend}". Choose from "min", "max", "both", "neither".')

        # Determine color ids for levels, i.e. position in 0-1 space
        # Length of these ids should be N + 1 -- that is, N - 1 colors
        # for values in-between levels, plus 2 colors for out-of-bounds.
        # * For same out-of-bounds colors, looks like [0, 0, ..., 1, 1]
        # * For unique out-of-bounds colors, looks like [0, X, ..., 1 - X, 1]
        #   where the offset X equals step/len(levels).
        # First get coordinates
        if not norm or type(norm) is mcolors.Normalize:
            norm = lambda x: x # linear transition, no need to normalize
        x_b = norm(levels)
        x_m = (x_b[1:] + x_b[:-1])/2 # get level centers after norm scaling
        y = (x_m - x_m.min())/(x_m.max() - x_m.min())
        if isinstance(y, ma.core.MaskedArray):
            y = y.filled(np.nan)
        y = y[np.isfinite(y)]
        # Account for out of bounds colors
        # WARNING: For some reason, "clipping" doesn't work when applying
        # norm, end up with unpredictable fill value if norm yields invalid
        # value. Must clip manually.
        offset = 0
        scale = 1
        eps = step/levels.size
        if extend in ('min','both'):
            offset = eps
            scale -= eps
        if extend in ('max','both'):
            scale -= eps
        y = np.concatenate(([0], offset + scale*y, [1])) # insert '0' (arg 3) before index '0' (arg 2)
        self._norm = norm
        self._x_b = x_b
        self._y = y
        if isinstance(norm, mcolors.LogNorm):
            self._norm_clip = (1e-250, None)
        else:
            self._norm_clip = None

        # Add builtin properties
        # NOTE: Are vmin/vmax even used?
        self.boundaries = levels
        self.vmin = levels.min()
        self.vmax = levels.max()
        self.clip = clip
        self.N = levels.size

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1."""
        # Follow example of LinearSegmentedNorm, but perform no interpolation,
        # just use searchsorted to bin the data.
        # Note the bins vector includes out-of-bounds negative (searchsorted
        # index 0) and out-of-bounds positive (searchsorted index N+1) values
        clip = self._norm_clip
        if clip:
            xq = np.clip(xq, *clip)
        xq = self._norm(np.atleast_1d(xq))
        yq = self._y[np.searchsorted(self._x_b, xq)] # which x-bin does each point in xq belong to?
        return ma.masked_array(yq, np.isnan(xq))

    def inverse(self, yq):
        """Raises error -- inversion after discretization is impossible."""
        raise RuntimeError('BinNorm is not invertible.')

#------------------------------------------------------------------------------#
# Normalizers intended to *pre-scale* levels passed to BinNorm
#------------------------------------------------------------------------------#
class LinearSegmentedNorm(mcolors.Normalize):
    """
    This is the default normalizer paired with `BinNorm` whenever `levels`
    are non-linearly spaced.

    It follows the example of the `~matplotlib.colors.LinearSegmentedColormap`
    source code and performs efficient, vectorized linear interpolation
    between the provided boundary levels. That is, the normalized value is
    linear with respect to its average **index** in the `levels` vector. This
    allows color transitions with uniform intensity across **arbitrarily
    spaced**, monotonically increasing points.

    Can be used by passing ``norm='segments'`` to any command accepting
    ``cmap``. The default midpoint is zero.
    """
    def __init__(self, levels, clip=False, **kwargs):
        """
        Parameters
        ----------
        levels : list of float
            The discrete data levels.
        **kwargs, clip
            Passed to `~matplotlib.colors.Normalize`.
        """
        # Save levels
        levels = np.atleast_1d(levels)
        if levels.size<=1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError(f'Levels {levels} passed to LinearSegmentedNorm must be monotonically increasing.')
        super().__init__(np.nanmin(levels), np.nanmax(levels), clip) # second level superclass
        self._x = levels
        self._y = np.linspace(0, 1, len(levels))

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1. Inverse operation
        of `~LinearSegmentedNorm.inverse`."""
        # Follow example of make_mapping_array for efficient, vectorized
        # linear interpolation across multiple segments.
        # Notes:
        # * Normal test puts values at a[i] if a[i-1] < v <= a[i]; for
        #   left-most data, satisfy a[0] <= v <= a[1]
        # * searchsorted gives where xq[i] must be inserted so it is larger
        #   than x[ind[i]-1] but smaller than x[ind[i]]
        x = self._x # from arbitrarily spaced monotonic levels
        y = self._y # to linear range 0-1
        xq = np.atleast_1d(xq)
        ind = np.searchsorted(x, xq)
        ind[ind==0] = 1
        ind[ind==len(x)] = len(x) - 1 # actually want to go to left of that
        distance = (xq - x[ind - 1])/(x[ind] - x[ind - 1])
        yq = distance*(y[ind] - y[ind - 1]) + y[ind - 1]
        return ma.masked_array(yq, np.isnan(xq))

    def inverse(self, yq):
        """Inverse operation of `~LinearSegmentedNorm.__call__`."""
        x = self._x
        y = self._y
        yq = np.atleast_1d(yq)
        ind = np.searchsorted(y, yq)
        ind[ind==0] = 1
        ind[ind==len(y)] = len(y) - 1
        distance = (yq - y[ind - 1])/(y[ind] - y[ind - 1])
        xq = distance*(x[ind] - x[ind - 1]) + x[ind - 1]
        return ma.masked_array(xq, np.isnan(yq))

class MidpointNorm(mcolors.Normalize):
    """
    Ensures a "midpoint" always lies at the central colormap color.
    Can be used by passing ``norm='midpoint'`` to any command accepting
    ``cmap``. The default midpoint is zero.

    Note
    ----
    See `this stackoverflow thread <https://stackoverflow.com/q/25500541/4970632>`_.
    """
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=None):
        """
        Parameters
        ----------
        midpoint : float, optional
            The midpoint, or the data value corresponding to the normalized
            value ``0.5`` -- halfway down the colormap.
        vmin, vmax, clip
            The minimum and maximum data values, and the clipping setting.
            Passed to `~matplotlib.colors.Normalize`.
        """
        # Bigger numbers are too one-sided
        super().__init__(vmin, vmax, clip)
        self._midpoint = midpoint

    def __call__(self, xq, clip=None):
        """Normalizes data values to the range 0-1. Inverse operation of
        `~MidpointNorm.inverse`."""
        # Get middle point in 0-1 coords, and value
        # Notes:
        # * Look up these three values in case vmin/vmax changed; this is
        #   a more general normalizer than the others. Others are 'parent'
        #   normalizers, meant to be static more or less.
        # * searchsorted gives where xq[i] must be inserted so it is larger
        #   than x[ind[i]-1] but smaller than x[ind[i]]
        #   x, y = [self.vmin, self._midpoint, self.vmax], [0, 0.5, 1]
        if self.vmin >= self._midpoint or self.vmax <= self._midpoint:
            raise ValueError(f'Midpoint {self._midpoint} outside of vmin {self.vmin} and vmax {self.vmax}.')
        x = np.array([self.vmin, self._midpoint, self.vmax])
        y = np.array([0, 0.5, 1])
        xq = np.atleast_1d(xq)
        ind = np.searchsorted(x, xq)
        ind[ind==0] = 1 # in this case will get normed value <0
        ind[ind==len(x)] = len(x) - 1 # in this case, will get normed value >0
        distance = (xq - x[ind - 1])/(x[ind] - x[ind - 1])
        yq = distance*(y[ind] - y[ind - 1]) + y[ind - 1]
        return ma.masked_array(yq, np.isnan(xq))
        # return ma.masked_array(np.interp(xq, x, y))

    def inverse(self, yq, clip=None):
        """Inverse operation of `~MidpointNorm.__call__`."""
        # Invert the above
        # x, y = [self.vmin, self._midpoint, self.vmax], [0, 0.5, 1]
        # return ma.masked_array(np.interp(yq, y, x))
        # Performs inverse operation of __call__
        x = np.array([self.vmin, self._midpoint, self.vmax])
        y = np.array([0, 0.5, 1])
        yq = np.atleast_1d(yq)
        ind = np.searchsorted(y, yq)
        ind[ind==0] = 1
        ind[ind==len(y)] = len(y) - 1
        distance = (yq - y[ind - 1])/(y[ind] - y[ind - 1])
        xq = distance*(x[ind] - x[ind - 1]) + x[ind - 1]
        return ma.masked_array(xq, np.isnan(yq))

#------------------------------------------------------------------------------#
# Register new colormaps; must come before registering the color cycles
# * If leave 'name' empty in register_cmap, name will be taken from the
#   Colormap instance. So do that.
# * Note that **calls to cmap instance do not interpolate values**; this is only
#   done by specifying levels in contourf call, specifying lut in get_cmap,
#   and using LinearSegmentedColormap.from_list with some N.
# * The cmap object itself only **picks colors closest to the "correct" one
#   in a "lookup table**; using lut in get_cmap interpolates lookup table.
#   See LinearSegmentedColormap doc: https://matplotlib.org/api/_as_gen/matplotlib.colors.LinearSegmentedColormap.html#matplotlib.colors.LinearSegmentedColormap
# * If you want to always disable interpolation, use ListedColormap. This type
#   of colormap instance will choose nearest-neighbors when using get_cmap, levels, etc.
#------------------------------------------------------------------------------#
def register_colors(nmax=np.inf, verbose=False):
    """
    Registers new color names and filters them to be
    sufficiently "perceptually distinct" in the HSL colorspace.
    Called on import.

    Use `~proplot.demos.color_show` to generate a table of the resulting
    filtered colors.
    """
    # Reset native colors dictionary and add some default groups
    # Add in CSS4 so no surprises for user, but we will not encourage this
    # usage and will omit CSS4 colors from the demo table.
    scale = (360, 100, 100)
    translate =  {'b': 'blue', 'g': 'green', 'r': 'red', 'c': 'cyan',
                  'm': 'magenta', 'y': 'yellow', 'k': 'black', 'w': 'white'}
    base = mcolors.BASE_COLORS
    full = {translate[key]:value for key,value in mcolors.BASE_COLORS.items()} # full names
    mcolors._colors_full_map.clear() # clean out!
    mcolors._colors_full_map.cache.clear() # clean out!
    for name,dict_ in (('base',base), ('full',full), ('css',mcolors.CSS4_COLORS)):
        colordict.update({name:dict_})

    # Register 'filtered' colors, and get their HSL values
    # The below order is order of preference for identical color names from
    # different groups.
    names = []
    seen = {*base, *full} # never overwrite these ones!
    hcls = np.empty((0,3))
    files = [os.path.join(_data_colors, f'{name}.txt') for name in ('opencolors', 'xkcd', 'crayola')]
    for file in files:
        category, _ = os.path.splitext(os.path.basename(file))
        data = np.genfromtxt(file, delimiter='\t', dtype=str, comments='%', usecols=(0,1)).tolist()
        # Immediately add all opencolors
        if category=='opencolors':
            dict_ = {name:color for name,color in data}
            mcolors._colors_full_map.update(dict_)
            colordict.update({'opencolors':dict_})
            continue
        # Other color dictionaries are filtered, and their names are sanitized
        i = 0
        dict_ = {}
        ihcls = []
        colordict[category] = {} # just initialize this one
        for name,color in data: # is list of name, color tuples
            if i>=nmax: # e.g. for xkcd colors
                break
            for regex,sub in _sanitize_names:
                name = re.sub(regex, sub, name)
            if name in seen or re.search(_bad_names, name):
                continue
            seen.add(name)
            names.append((category, name)) # save the category name pair
            ihcls.append(to_xyz(color, space=_distinct_colors_space))
            dict_[name] = color # save the color
            i += 1
        _colors_unfiltered[category] = dict_
        hcls = np.concatenate((hcls, ihcls), axis=0)

    # Remove colors that are 'too similar' by rounding to the nearest n units
    # WARNING: Unique axis argument requires numpy version >=1.13
    deleted = 0
    hcls = hcls/np.array(scale)
    hcls = np.round(hcls/_distinct_colors_threshold).astype(np.int64)
    _, index, counts = np.unique(hcls, return_index=True, return_counts=True, axis=0) # get unique rows
    counts = counts.sum()
    for i,(category,name) in enumerate(names):
        if name not in _exceptions_names and i not in index:
            deleted += 1
        else:
            colordict[category][name] = _colors_unfiltered[category][name]
    for key,kw in colordict.items():
        mcolors._colors_full_map.update(kw)
    if verbose:
        print(f'Started with {len(names)} colors, removed {deleted} insufficiently distinct colors.')

def register_cmaps():
    """
    Registers colormaps packaged with ProPlot or added by the user. That is,
    add maps to the `matplotlib.cm.cmap_d` dictionary. Called on import.

    Use `~proplot.demos.cmap_show` to generate a table of the resulting
    color cycles.
    """
    # First read from file
    N_hires = rcParams['image.lut']
    for filename in sorted(glob.glob(os.path.join(_data_cmaps, '*'))) + \
            sorted(glob.glob(os.path.join(_data_user, '*'))):
        # Read table of RGB values
        if not re.search('\.(x?rgba?|json|xml)$', filename):
            continue
        name = os.path.basename(filename)
        name = name.split('.')[0]
        # if name in mcm.cmap_d:
        #     continue
        # Read .rgb, .rgba, .xrgb, and .xrgba files
        if re.search('\.x?rgba?$', filename):
            # Load
            ext = filename.split('.')[-1]
            try:
                cmap = np.loadtxt(filename, delimiter=',') # simple
            except:
                print(f'Failed to load {os.path.basename(filename)}.')
                continue
            # Build x-coordinates and standardize shape
            N = cmap.shape[0]
            if ext[0] != 'x':
                x = np.linspace(0, 1, N)
                cmap = np.concatenate((x[:,None], cmap), axis=1)
            if cmap.shape[1] not in (4,5):
                raise ValueError(f'Invalid number of columns for colormap "{name}": {cmap.shape[1]}.')
            if (cmap[:,1:4]>10).any(): # from 0-255 to 0-1
                cmap[:,1:4] = cmap[:,1:4]/255
            # Build color dict
            x = cmap[:,0]
            x = (x - x.min()) / (x.max() - x.min()) # for some reason, some aren't in 0-1 range
            if cmap.shape[1]==5:
                channels = ('red', 'green', 'blue', 'alpha')
            else:
                channels = ('red', 'green', 'blue')
            # Optional cycles
            if re.match('(cycle|qual)[0-9]', name.lower()):
                cmap = mcolors.ListedColormap(cmap[:,1:])
                cycles.add(name)
            else:
                cdict = {}
                for i,channel in enumerate(channels):
                    vector = cmap[:,i+1:i+2]
                    cdict[channel] = np.concatenate((x[:,None], vector, vector), axis=1).tolist()
                cmap = mcolors.LinearSegmentedColormap(name, cdict, N) # using static method is way easier
                cmaps.add(name)
        # Load XML files created with scivizcolor
        # Adapted from script found here: https://sciviscolor.org/matlab-matplotlib-pv44/
        elif re.search('\.xml$', filename):
            try:
                xmldoc = etree.parse(filename)
            except IOError:
                raise ValueError('The input file is invalid. It must be a colormap xml file. Go to https://sciviscolor.org/home/colormaps/ for some good options.')
            x = []
            colors = []
            for s in xmldoc.getroot().findall('.//Point'):
                x.append(float(s.attrib['x']))
                colors.append((float(s.attrib['r']), float(s.attrib['g']), float(s.attrib['b'])))
            N = len(x)
            x = np.array(x)
            x = (x - x.min()) / (x.max() - x.min()) # for some reason, some aren't in 0-1 range
            colors = np.array(colors)
            if re.match('(cycle|qual)[0-9]?', name.lower()):
                cmap = mcolors.ListedColormap([to_rgb(color) for color in colors])
                cycles.add(name)
            else:
                cdict = {}
                for i,channel in enumerate(('red', 'green', 'blue')):
                    vector = colors[:,i:i+1]
                    cdict[channel] = np.concatenate((x[:,None], vector, vector), axis=1).tolist()
                cmap = mcolors.LinearSegmentedColormap(name, cdict, N) # using static method is way easier
                cmaps.add(name)
        # Directly read segmentdata of hex strings
        # Will ensure that HSL colormaps have the 'space' entry
        else:
            with open(filename, 'r') as file:
                segmentdata = json.load(file)
            if 'space' in segmentdata:
                space = segmentdata.pop('space')
                cmap = PerceptuallyUniformColormap(name, segmentdata, space=space, N=N_hires)
            else:
                cmap = mcolors.LinearSegmentedColormap(name, segmentdata, N=N_hires)
            cmaps.add(name)
        # Register maps (this is just what register_cmap does)
        # If the _r (reversed) version is stored on file, store the straightened one
        if re.search('_r$', name):
            name = name[:-2]
            cmap = cmap.reversed()
            cmap.name = name
        mcm.cmap_d[name] = cmap

    # Fix the builtin rainbow colormaps by switching from Listed to
    # LinearSegmented -- don't know why matplotlib shifts with these as
    # discrete maps by default, dumb.
    for name in _cmap_categories['Matplotlib Originals']: # initialize as empty lists
        cmap = mcm.cmap_d.get(name, None)
        if cmap and isinstance(cmap, mcolors.ListedColormap):
            mcm.cmap_d[name] = mcolors.LinearSegmentedColormap.from_list(name, cmap.colors)

    # Reverse some included colormaps, so colors
    # go from 'cold' to 'hot'
    for name in ('Spectral',):
        mcm.cmap_d[name] = mcm.cmap_d[name].reversed()

    # Add shifted versions of cyclic colormaps, and prevent same colors on ends
    # TODO: Add automatic shifting of colormap by N degrees in the CmapDict
    for name in ['twilight', 'Phase', 'GrayCycle']:
        cmap = mcm.cmap_d.get(name, None)
        if not isinstance(cmap, mcolors.LinearSegmentedColormap):
            continue
        # Shift data
        data = cmap._segmentdata
        data_shift = data.copy()
        for key,array in data.items():
            array = np.array(array)
            # Drop an end color
            array = array[1:,:]
            array_shift = array.copy()
            array_shift[:,0] -= 0.5
            array_shift[:,0] %= 1
            array_shift = array_shift[array_shift[:,0].argsort(),:]
            # Normalize x-range
            array[:,0] -= array[:,0].min()
            array[:,0] /= array[:,0].max()
            data[key] = array
            array_shift[:,0] -= array_shift[:,0].min()
            array_shift[:,0] /= array_shift[:,0].max()
            data_shift[key] = array_shift
        # Register shifted version and original
        for i,data in enumerate((data, data_shift)):
            name = name if i==0 else f'{name}_shifted'
            cmap = mcolors.LinearSegmentedColormap(name, data, cmap.N)
            cmap._cyclic = True
            mcm.cmap_d[name] = cmap

    # Delete ugly cmaps (strong-arm user into using the better ones)
    # TODO: Better way to generalize this language stuff? Not worth it maybe.
    greys = mcm.cmap_d.get('Greys', None)
    if greys is not None:
        mcm.cmap_d['Grays'] = greys
    # TODO: Add this to cmap dict __init__?
    for category in _cmap_categories_delete:
        for name in _cmap_categories:
            mcm.cmap_d.pop(name, None)

def register_cycles():
    """
    Registers color cycles defined by ``.hex`` files (each of which contains a
    comma-delimited list of HEX strings) packaged with ProPlot or
    added by the user. Also registers some cycles hardcoded into this module.
    Called on import.

    Use `~proplot.demos.cycle_show` to generate a table of the resulting
    color cycles.
    """
    # Read lists of hex strings from disk
    for filename in sorted(glob.glob(os.path.join(_data_cmaps, '*'))) + \
            sorted(glob.glob(os.path.join(_data_user, '*'))):
        if not re.search('\.hex$', filename):
            continue
        name = os.path.basename(filename)
        name = name.split('.hex')[0]
        colors = ''.join(s.strip() for s in open(filename)).split(',') # single or multiple lines
        if len(colors)<2:
            raise ValueError(f'Error reading file "{filename}".')
        colors = [mcolors.to_rgb(c) for c in colors] # from list of tuples
        _cycles_loaded[name] = colors

    # Register names
    # Note that 'loaded' cycles will overwrite any presets with same name
    for name,colors in {**_cycles_preset, **_cycles_loaded}.items():
        mcm.cmap_d[name] = mcolors.ListedColormap([to_rgb(color) for color in colors], name=name)
        cycles.add(name)

    # Remove some redundant or ugly ones
    for key in ('tab10', 'tab20', 'Paired', 'Pastel1', 'Pastel2', 'Dark2'):
        mcm.cmap_d.pop(key, None)
    # *Change* the name of some more useful ones
    for (name1,name2) in [('Accent','Set1'), ('tab20b','Set4'), ('tab20c','Set5')]:
        orig = mcm.cmap_d.pop(name1, None)
        if orig:
            mcm.cmap_d[name2] = orig
            cycles.add(name2)

# Register stuff when this module is imported
# The 'cycles' are simply listed colormaps, and the 'cmaps' are the smoothly
# varying LinearSegmentedColormap instances or subclasses thereof
cmaps = set() # track *downloaded* colormaps; user can then check this list
"""List of new registered colormap names."""

cycles = set() # track *all* color cycles
"""List of registered color cycle names."""

_colors_unfiltered = {} # downloaded colors categorized by filename
colordict = {} # limit to 'sufficiently unique' color names
"""Filtered, registered color names by category."""

register_colors() # must be done first, so we can register OpenColor cmaps
register_cmaps()
register_cycles()
cmaps = set(sorted(cmaps))
cycles = set(sorted(cycles))

# Finally our dictionary of normalizers
# Includes some custom classes, so has to go at end
# NOTE: Make BinNorm inaccessible to users. Idea is that all other normalizers
# can be wrapped by BinNorm -- BinNorm is just used to break colors into
# discrete levels.
normalizers = {
    'none':       mcolors.NoNorm,
    'null':       mcolors.NoNorm,
    'zero':       MidpointNorm,
    'midpoint':   MidpointNorm,
    'segments':   LinearSegmentedNorm,
    'segmented':  LinearSegmentedNorm,
    'log':        mcolors.LogNorm,
    'linear':     mcolors.Normalize,
    'power':      mcolors.PowerNorm,
    'symlog':     mcolors.SymLogNorm,
    }
"""Dictionary of possible normalizers. See `Norm` for a table."""

