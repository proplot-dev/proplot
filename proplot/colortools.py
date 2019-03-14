#!/usr/bin/env python3
"""
Registers colormaps, color cycles, and color string names with
`register_cmaps`, `register_cycles`, and `register_colors`.
Defines tools for creating new colormaps and color cycles, i.e. `colormap`
and `colors`. Defines helpful new `~matplotlib.colors.Normalizer` and
`~matplotlib.colors.Colormap` classes.

.. raw:: html

   <h1>Perceptually uniform colormaps</h1>


ProPlot’s custom colormaps, and its `PerceptuallyUniformColormap` class,
consist of *linear transitions* between channel values in any of
three possible "perceptually uniform", HSV-like colorspaces.
These colorspaces can be described as follows (see also `these 
visualizations <http://www.hsluv.org/comparison/>`_):

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

.. raw:: html

   <h1>From other projects</h1>


I’ve removed some outdated “miscellaneous” colormaps that are packaged
by default (see `this reference
<https://matplotlib.org/examples/color/colormaps_reference.html>`_),
added some custom colormaps generated with the `PerceptuallyUniformColormap`
class, and added new "perceptually uniform" colormaps from the following projects:

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


.. raw:: html

   <h1>Flexible colormap declaration</h1>


All of the `~matplotlib.axes.Axes` methods listed in `cmap_methods` and
`cycle_methods` have been wrapped. For the latter methods, a brand new
keyword arg called `cycle` has been added, for changing the axes property
cycle on-the-fly. The `cmap` and `cycle` arguments are all passed through
the magical `Colormap` function.

`Colormap` is incredibly powerful -- it can make colormaps on-the-fly, look
up existing maps, and merge them. As such, any of the following are now
valid `cmap` and `cycle` arguments:

1. Registered colormap names. For example, ``'Blues'`` or ``'Sunset'``.
   See :ref:`Colors` for a visalization of the registered maps.

   Cycles are constructed automatically be sampling colors from these colormaps;
   use e.g. ``('Blues', 5)`` to specify the number of colors in the cycle.
2. Gradations of a single hue. For example, ``'maroon'`` creates colormap spanning
   white to maroon, and ``'maroon90'`` spans from a 90% luminance pale red
   color to maroon (see `~proplot.colors.Colormap`).

   Cycles are once again constructed automatically; use e.g. ``('maroon', 10)``
   to specify the number of colors in the cycle.
3. Registered color cycle names. For example, ``'Cycle1'``. See :ref:`Colors`
   for a visualization of the registered cycles.
4. Lists of colors. For example, ``['red', 'blue', 'green']``.
   Cycles are constructed automatically from these
5. Dictionary containing the keys ``'h'``, ``'s'``, and ``'l'``. This builds
   a `PerceptuallyUniformColormap` using the
   `~PerceptuallyUniformColormap.from_hsl` constructor.
6. **List of any of the above five arguments, to merge the resulting
   colormaps.** For
   example, use ``['viridis', 'navy']`` to merge the virids map with a colormap
   spanning navy to white.

Note when assigning to the ``proplot.rc.cycle`` global setting (see
`~proplot.rcmod`), the argument is also interpreted as above. For example,
``proplot.rc.cycle = ('blue', 10)`` will construct a color cycle with 10 colors
ranging from white to blue.
"""
#------------------------------------------------------------------------------#
# Notes
#------------------------------------------------------------------------------#
# Note colormaps are *callable*, will just deliver the corresponding color, easy.
# Notes on different colorspaces:
# TODO: Allow some colormaps (e.g. topography) to have ***fixed*** levels,
# i.e. force the user to use very high resolution (e.g. 256 levels) and the
# 'x' coordinates will enforce discrete jumps between colors. Consider doing
# this for 'seismic' map, and others.
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
from matplotlib import rcParams
from . import utils, colormath
from .utils import _default, ic
_data_user = os.path.join(os.path.expanduser('~'), '.proplot')
_data_cmaps = os.path.join(os.path.dirname(__file__), 'cmaps') # or parent, but that makes pip install distribution hard
_data_colors = os.path.join(os.path.dirname(__file__), 'colors') # or parent, but that makes pip install distribution hard

# Default number of colors
_N_hires = 256

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
# _distinct_colors_threshold = 0.07
_distinct_colors_threshold = 0.07 # bigger number equals fewer colors
_distinct_colors_space = 'hsl' # register colors distinct in this space?
_distinct_colors_exceptions = [
    'white', 'black', 'gray', 'red', 'pink', 'grape',
    'sky blue', 'eggshell', 'sea blue',
    'violet', 'indigo', 'blue',
    'coral', 'tomato red', 'crimson',
    'cyan', 'teal', 'green', 'lime', 'yellow', 'orange',
    'red orange', 'yellow orange', 'yellow green', 'blue green',
    'blue violet', 'red violet',
    ]
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

# Names of builtin colormaps
# NOTE: Has support for 'x' coordinates in first column.
# NOTE: For 'alpha' column, must use a .rgba filename
# TODO: Better way to save colormap files.
_cmap_categories = { # initialize as empty lists
    # We keep these ones
    'Matplotlib Originals': [
        'viridis', 'plasma', 'inferno', 'magma', 'twilight', 'twilight_shifted',
        ],

    # Assorted origin, but these belong together
    'Grayscale': [
        'Grays',
        'GrayCM',
        'GrayC',
        'PseudoGray',
        'GrayCycle',
        'GrayCycle_shifted',
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
        # 'GYPi', 'GnRP', 'BrBG', 'PuOr', 'GyRd', 'BuRd', 'BuYlRd', 'GnYlRd', 'Spectral'
        'Spectral', 'Spectral_cut', 'PiYG', 'PiYG_cut', 'PRGn', 'PRGn_cut',
        'BrBG', 'BrBG_cut', 'PuOr', 'PuOr_cut', 'RdGY', 'RdGY_cut',
        'RdBu', 'RdBu_cut', 'RdYlBu', 'RdYlBu_cut', 'RdYlGn', 'RdYlGn_cut',
        ],

    # Custom maps
    'ProPlot Sequential': [
         'Glacial',
        'Bog', 'Verdant',
        'Lake', 'Turquoise', 'Forest',
        'Blood',
        'Sunrise', 'Sunset', 'Fire',
        'Golden'
        ],
        # 'Vibrant'], # empty at first, fill automatically
    'ProPlot Diverging': [
        'IceFire', 'NegPos', 'BlueRed', 'PurplePink', 'DryWet', 'DrierWetter', 'LandSea'
        ],

    # Other
    # BlackBody2 is actually from Los Alamos, and a couple are from Kenneth's
    # website, but organization is better this way.
    'Miscellaneous': [
        'Temperature', # from ???
        'BWR',
        'ColdHot',
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

    # cmOcean
    'cmOcean Sequential': [
        'Oxy', 'Thermal', 'Dense', 'Ice', 'Haline',
        'Deep', 'Algae', 'Tempo', 'Speed', 'Turbid', 'Solar', 'Matter',
        'Amp', 'Phase', 'Phase_shifted'
        ],
    'cmOcean Diverging': [
        'Balance', 'Curl', 'Delta'
        ],

    # Statistik
    'Statistik Stadt Zürich': [
        'MutedBlue', 'MutedRed', 'MutedDry', 'MutedWet',
        'MutedBuRd', 'MutedBuRd_cut', 'MutedDryWet', 'MutedDryWet_cut',
        ],

    # Kenneth Moreland
    # See: http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/index.html
    # Soft coolwarm from: https://www.kennethmoreland.com/color-advice/
    # 'Kenneth Moreland Sequential': [
    #     'BlackBody', 'Kindlmann', 'ExtendedKindlmann',
    #     ],
    'Kenneth Moreland': [
        'CoolWarm', 'MutedCoolWarm', 'SoftCoolWarm',
        'BlueTan', 'PurpleOrange', 'CyanMauve', 'BlueYellow', 'GreenRed',
        ],

    # CET isoluminant maps
    # See: https://peterkovesi.com/projects/colourmaps/
    # All the others have better options
    'Isoluminant': [
        'Iso1', 'Iso2', 'Iso3', 
        # 'CET1', 'CET2', 'CET3', 'CET4',
        ],
    # 'CET Rainbow': [
    #     ],
    # 'CET Diverging': [
    #     ],
    # 'CET Cyclic': [
    #     ],

    # Sky themes from rafi; not very scientifically useful, but pretty
    # 'Sky' : [
    #     'Sky1', 'Sky2', 'Sky3', 'Sky4', 'Sky5', 'Sky6', 'Sky7',
    #     ],

    # Los Alamos
    # See: https://datascience.lanl.gov/colormaps.html
    # Most of these have analogues in SciVisColor, added the few unique
    # ones to Miscellaneous category
    # 'Los Alamos Sequential': [
    #     'MutedRainbow', 'DarkRainbow', 'MutedBlue', 'DeepBlue', 'BrightBlue', 'BrightGreen', 'WarmGray',
    #     ],
    # 'Los Alamos Diverging': [
    #     'MutedBlueGreen', 'DeepBlueGreen', 'DeepBlueGreenAsym', 'DeepColdHot', 'DeepColdHotAsym', 'ExtendedCoolWarm'
    #     ],

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

    # Remove the "combo" maps (and ugly diverging ones) because these can
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

    }


# Categories to ignore/*delete* from dictionary because they suck donkey balls
_cmap_categories_delete = ['Alt Diverging', 'Alt Sequential', 'Alt Rainbow', 'Miscellaneous Orig']

# Slice indices that split up segments of names
# WARNING: Must add to this list manually! Not worth trying to generalize.
# List of string cmap names, and the indices where they can be broken into parts
_cmap_parts = {
    # Sequential
    # Decided these shouldn't be reversed; left colors always refer
    # to 'light' colors, right to 'dark' colors. Also there is BuPu and PuBu
    # 'ylorbr':       (None, 2, 4, None),
    # 'ylorrd':       (None, 2, 4, None),
    # 'orrd':         (None, 2, None),
    # 'purd':         (None, 2, None),
    # 'rdpu':         (None, 2, None),
    # 'bupu':         (None, 2, None),
    # 'gnbu':         (None, 2, None),
    # 'pubu':         (None, 2, None),
    # 'ylgnbu':       (None, 2, 4, None),
    # 'pubugn':       (None, 2, 4, None),
    # 'bugn':         (None, 2, None),
    # 'ylgn':         (None, 2, None),
    # Diverging
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
# Special class for colormap names
#------------------------------------------------------------------------------#
class CmapDict(dict):
    """
    Flexible, case-insensitive colormap identification. Behaves like a
    dictionary, with three new features:

    1. Names are case insensitive: ``'Blues'``, ``'blues'``, and ``'bLuEs'``
       are all valid names for the "Blues" colormap.
    2. "Reversed" colormaps are not stored directly: Requesting e.g.
       ``'Blues_r'`` will just look up ``'Blues'``, then return the result
       of the `~matplotlib.colors.Colormap.reversed` method.
    3. Diverging colormap names can be referenced by their "inverse" name.
       For example, ``'BuRd'`` is equivalent to ``'RdBu_r'``, as are
       ``'BuYlRd'`` and ``'RdYlBu_r'``.

    ProPlot replaces the `matplotlib.cm.cmap_d` dictionary with an instance
    of `CmapDict`.
    """
    # Initialize -- converts keys to lower case and
    # ignores the 'reverse' maps
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
        # Call this to avoid sanitization
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
        # Assume lowercase
        key = self._sanitize_key(key)
        return self._getitem(key)

    def __setitem__(self, key, item):
        # Set item
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key}. Must be string.')
        return super().__setitem__(key.lower(), item)

    def __contains__(self, item):
        # Must be overriden?
        try:
            self.__getitem__(item)
            return True
        except KeyError:
            return False

    # Other methods
    def get(self, key, *args):
        """
        Case-insensitive version of ``dict.get``.
        """
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
        """
        Case-insensitive version of ``dict.pop``.
        """
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

# Override entire colormap dictionary
if not isinstance(mcm.cmap_d, CmapDict):
    mcm.cmap_d = CmapDict(mcm.cmap_d)

#------------------------------------------------------------------------------#
# More generalized utility for retrieving colors
#------------------------------------------------------------------------------#
def _get_space(space):
    """
    Verify requested colorspace is valid.
    """
    space = _space_aliases.get(space, None)
    if space is None:
        raise ValueError(f'Unknown colorspace "{space}".')
    return space

def shade(color, shade=0.5):
    """
    Change the "shade" of a color by messing with its luminance channel.
    """
    try:
        color = mcolors.to_rgb(color) # ensure is valid color
    except Exception:
        raise ValueError(f'Invalid RGBA argument {color}. Registered colors are: {", ".join(mcolors._colors_full_map.keys())}.')
    color = [*colormath.rgb_to_hsl(*color)]
    color[2] = max([0, min([color[2]*shade, 100])]) # multiply luminance by this value
    color = [*colormath.hsl_to_rgb(*color)]
    return tuple(color)

def to_rgb(color, space='rgb'):
    """
    Generalization of `~matplotlib.colors.to_rgb`. Translates color tuples
    from *any* colorspace to rgb. Also will convert color strings to tuple.
    """
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
    elif space=='rgb':
        color = color[:3] # trim alpha
        if any(c>1 for c in color):
            color = [c/255 for c in color] # scale to within 0-1
    else:
        raise ValueError('Invalid RGB value.')
    return color

def to_xyz(color, space):
    """
    Inverse of above, translate to some colorspace.
    """
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

def add_alpha(color):
    """
    Ensures presence of alpha channel.
    """
    if not np.iterable(color) or isinstance(color, str):
        raise ValueError('Input must be color tuple.')
    if len(color)==3:
        color = [*color, 1.0]
    elif len(color)==4:
        color = [*color] # copy, and put into list
    else:
        raise ValueError(f'Tuple length must be 3 or 4, got {len(color)}.')
    return color

def get_channel_value(color, channel, space='hsl'):
    """
    Gets hue, saturation, or luminance channel value from registered
    string color name.

    Arguments
    ---------
        color :
            scalar numeric ranging from 0-1, or string color name, optionally
            with offset specified as '+x' or '-x' at the end of the string for
            arbitrary float x.
        channel :
            channel number or name (e.g., 0, 1, 2, 'h', 's', 'l')
    """
    # Interpret channel
    channel_idxs = {'hue': 0, 'saturation': 1, 'chroma': 1, 'luminance': 2,
                    'alpha': 3, 'h': 0, 's': 1, 'c': 1, 'l': 2}
    channel = channel_idxs.get(channel, channel)
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

#------------------------------------------------------------------------------#
# Generalized colormap/cycle constructors
#------------------------------------------------------------------------------#
def clip_cmap(cmap, left=None, right=None, N=None):
    """
    Helper function that cleanly divides linear segmented colormaps and
    subsamples listed colormaps.

    Parameters
    ----------
    cmap : `~matplotlib.colors.Colormap`
        The colormap.
    left : None or float, optional
        Clips the "left" of the colormap. For example, ``left=0.1``
        deletes the leftmost 10%.

        For `~matplotlib.colors.ListedColormap` maps, this is an integer index
        in the color list.
    right : None or float, optional
        Clips the "right" of the colormap. For example, ``right=0.9``
        deletes the rightmost 10%.

        For `~matplotlib.colors.ListedColormap` maps, this is an integer index
        in the color list.
    N : None or int, optional
        For `~matplotlib.colors.ListedColormap` maps, this is an alias for
        `right`.  Otherwise, it is ignored.
    """
    # Optionally clip edges or resample map.
    if isinstance(cmap, mcolors.ListedColormap):
        slicer = None
        if left is not None or right is not None:
            slicer = slice(left, right)
        elif N is not None:
            slicer = slice(None, N)
        # Just sample indices for listed maps
        if slicer:
            try:
                cmap = mcolors.ListedColormap(cmap.colors[slicer])
            except Exception:
                raise ValueError(f'Invalid indices {slicer} for listed colormap.')
    elif left is not None or right is not None:
        # Trickier for segment data maps
        # First get segmentdata and parse input
        kwargs = {}
        olddata = cmap._segmentdata
        newdata = {}
        if left is None:
            left = 0
        if right is None:
            right = 1
        if hasattr(cmap, 'space'):
            kwargs['space'] = cmap.space
        # Next resample the segmentdata arrays
        for key,xyy in olddata.items():
            if key in ('gamma1', 'gamma2'):
                newdata[key] = xyy
                continue
            xyy = np.array(xyy)
            x = xyy[:,0]
            xleft, = np.where(x>left)
            xright, = np.where(x<right)
            if len(xleft)==0:
                raise ValueError(f'Invalid x minimum {left}.')
            if len(xright)==0:
                raise ValueError(f'Invalid x maximum {right}.')
            l, r = xleft[0], xright[-1]
            newxyy = xyy[l:r+1,:].copy()
            if l>0:
                xl = xyy[l-1,1:] + (left - x[l-1])*(xyy[l,1:] - xyy[l-1,1:])/(x[l] - x[l-1])
                newxyy = np.concatenate(([[left, *xl]], newxyy), axis=0)
            if r<len(x)-1:
                xr = xyy[r,1:] + (right - x[r])*(xyy[r+1,1:] - xyy[r,1:])/(x[r+1] - x[r])
                newxyy = np.concatenate((newxyy, [[right, *xr]]), axis=0)
            newxyy[:,0] = (newxyy[:,0] - left)/(right - left)
            newdata[key] = newxyy
        # And finally rebuild map
        cmap = type(cmap)(cmap.name, newdata, **kwargs)
    return cmap

def Colormap(*args, name=None, N=None, 
        extend='both',
        left=None, right=None, x=None, reverse=False, # optionally truncate color range by these indices
        ratios=1, gamma=None, gamma1=None, gamma2=None,
        register=True, save=False,
        **kwargs):
    """
    Convenience function for generating and **merging** colormaps
    in a variety of ways.

    Parameters
    ----------
    *args : `~matplotib.colors.Colormap`, dict, list of str, or str
        Each arg generates a single colormap. If ``len(args)>1``, the colormaps
        are merged.

        If the arg is a `~matplotlib.colors.Colormap`, nothing more is done.
        Otherwise, the colormap is generated as follows:

        * If arg is dict, it is passed to `~PerceptuallyUniformColormap.from_hsl`.
        * If arg is a list of str, it is assumed a list of color names or
          hex names, and is used to make a `~matplotlib.colors.ListedColormap`.
        * If arg is a str and is a "registered" colormap name (i.e. is in
          the `matplotlib.cm.cmap_d` dictionary), that colormap is used.
        * If arg is a str and is a "registered" color name (i.e. is in
          the `matplotlib.colors._colors_full_map` dictionary), a
          monochromatic colormap is generated with `monochrome_cmap`, spanning
          from 100% luminance and 0% saturation to that color.

          The color string can also look like ``'name90'``, where the trailing
          number indicates the luminance at the other end of the colormap (by
          default, this is 100).

    name : None or str, optional
        Name of colormap. Default name is ``'no_name'``.

        The resulting colormap can then be invoked by passing ``cmap='name'``
        to plotting functions like `~matplotlib.figure.Figure.contourf`.
    N : None or int, optional
        Number of colors to generate in the hidden lookupt table ``_lut``.
        By default, a relatively high resolution of 256 is chosen (see notes).
    extend : {'both', 'min', 'max', 'neither'}, optional
        Specifies sides for which you want "out-of-bounds" data to have
        their own color. This improves upon the matplotlib API by ensuring
        **the colors on either end of the colorbar are always the most intense
        colors in the colormap** -- no matter the extend property. By default
        when `extend` is not ``'both'``, matplotlib just lobs off the most
        intense colors on either end.

        Note you can still use ``extend='neither'`` in your call to `Colormap`
        with ``extend='both'`` in the contour or colorbar call. This
        means that colors at the ends of the main region will be same as
        the out-of-bounds colors.
    left, right : float or list of float
        Optionally *delete* colors on the left and right sides of the
        colormap(s). For example, ``left=0.1`` deletes the leftmost
        10% of the colormap; ``right=0.9`` deletes the rightmost 10%.

        If list, length must match ``len(args)``, and applies to *each*
        colormap in the list before they are merged. If float, applies
        to the *final* colormap. No difference if ``len(args)`` is 1.
    x : (float, float), optional
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

        So far no way to apply this differently for each colormap in `args`.
    save : bool, optional
        Whether to save the colormap in the folder ``~/.proplot``. The
        folder is created if it does not already exist.

        If the colormap is a `~matplotlib.colors.ListedColormap` (i.e. a
        "color cycle"), the list of hex strings are written to ``name.hex``.

        If the colormap is a `~matplotlib.colors.LinearSegmentedColormap`,
        the segment data dictionary is written to ``name.json``.

    Notes
    -----
    Resampling method described in `this post
    <https://stackoverflow.com/q/48613920/4970632>`_ is awful! All it does is
    reduce the lookup table size -- what ends up happening under the hood is
    matplotlib tries to *evenly* draw ``N-1`` (``'min'`` or ``'max'``) or ``N-2``
    (``'neither'``) colors from a lookup table with ``N`` colors, which means
    it simply *skips over* 1 or 2 colors in the middle of the lookup table,
    which will cause visual jumps!

    Instead, we make stored lookup table completely independent from the
    number of levels; will have many high-res segments while the colormap
    `N` is very small.

    Warning
    -------
    Does colormap saving work for non-perceptually uniform maps or
    combinations thereof?
    """
    # Turns out pcolormesh makes QuadMesh, which itself is a Collection,
    # which itself gets colors when calling draw() using update_scalarmappable(),
    # which itself uses to_rgba() to get facecolors, which itself is an inherited
    # ScalarMappable method that simply calls the colormap with numbers.
    # Anyway the issue *has* to be with pcolor, because when giving pcolor an
    # actual instance, no longer does that thing where final levels equal extensions.
    # Since collection API does nothing to underlying data or cmap, must be
    # something done by pcolormesh function.
    _N = N or _N_hires
    imaps = []
    name = name or 'no_name' # must have name, mcolors utilities expect this
    if len(args)==0:
        args = [rcParams['image.cmap']] # use default

    for i,cmap in enumerate(args):
        # Retrieve Colormap instance
        # Also make sure you reset the lookup table (get_cmap does this
        # by calling _resample).
        # TODO: What if colormap names conflict with color names! Maybe
        # address this! Currently, this makes it impossible to make a
        # monochrome colormap from some named color if that name also
        # exists for a colormap.
        if cmap is None:
            cmap = rcParams['image.cmap']
        if isinstance(cmap,str) and cmap in mcm.cmap_d:
            cmap = mcm.cmap_d[cmap]
            if isinstance(cmap, mcolors.LinearSegmentedColormap):
                cmap = cmap._resample(_N)
        if isinstance(cmap, mcolors.Colormap):
            # Allow overriding the gamma, otherwise do nothing
            if isinstance(cmap, PerceptuallyUniformColormap):
                if gamma and not gamma1 and not gamma2:
                    gamma1 = gamma2 = gamma
                if gamma1 or gamma2:
                    segmentdata = cmap._segmentdata.copy()
                    if gamma1:
                        segmentdata['gamma1'] = gamma1
                    if gamma2:
                        segmentdata['gamma2'] = gamma2
                    cmap = type(cmap)(cmap.name, segmentdata, space=cmap.space, mask=cmap.mask)
            elif isinstance(cmap, mcolors.LinearSegmentedColormap):
                if gamma:
                    cmap._gamma = gamma
                    cmap._init()
        elif isinstance(cmap, dict):
            # Dictionary of hue/sat/luminance values or 2-tuples representing linear transition
            save = cmap.pop('save', save)
            name = cmap.pop('name', name)
            for key in cmap:
                if key in kwargs:
                    warnings.warn(f'Got duplicate keys "{key}" in cmap dictionary ({cmap[key]}) and in keyword args ({kwargs[key]}). Using first one.')
            kw = kwargs.update
            cmap = PerceptuallyUniformColormap.from_hsl(name, N=_N, **{**kwargs, **cmap})
        elif not isinstance(cmap, str):
            # List of colors
            cmap = mcolors.ListedColormap(cmap, name=name, **kwargs)
        else:
            # Monochrome colormap based from input color (i.e. single hue)
            regex = '([0-9].)$'
            match = re.search(regex, cmap) # declare maximum luminance with e.g. red90, blue70, etc.
            cmap = re.sub(regex, '', cmap) # remove options
            fade = kwargs.pop('fade', 90) if not match else float(match.group(1)) # default fade to 100 luminance
            # Build colormap
            cmap = to_rgb(cmap) # to ensure is hex code/registered color
            cmap = monochrome_cmap(cmap, fade, name=name, N=_N, **kwargs)

        # Transform each colormap
        left_i, right_i = None, None
        if np.iterable(left):
            left_i = left[i]
        if np.iterable(right):
            right_i = right[i]
        cmap = clip_cmap(cmap, left_i, right_i, N=N)
        if np.iterable(reverse) and reverse[i]:
            cmap = cmap.reversed()

        # Add to list
        imaps += [cmap]

    # Now merge the result of this arbitrary user input
    # Since we are merging cmaps, potentially *many* color transitions; use big number by default
    if len(imaps)>1:
        N_merge = _N*len(imaps)
        cmap = merge_cmaps(*imaps, name=name, ratios=ratios, N=N_merge)

    # Transform merged colormap
    if np.iterable(left): # was applied to each map, not the final one
        left = None
    if np.iterable(right):
        right = None
    cmap = clip_cmap(cmap, left, right, N=N)
    if not np.iterable(reverse) and reverse:
        cmap = cmap.reversed()

    if isinstance(cmap, mcolors.LinearSegmentedColormap) and N is not None:
        # Perform a crude resampling of the data, i.e. just generate a
        # low-resolution lookup table instead.
        # NOTE: All this does is create a new colormap with *attribute* N levels,
        # for which '_lut' attribute has not been generated yet.
        offset = {'neither':-1, 'max':0, 'min':0, 'both':1}
        if extend not in offset:
            raise ValueError(f'Unknown extend option {extend}.')
        cmap = cmap._resample(N - offset[extend]) # see mcm.get_cmap source

    # Register the colormap
    mcm.cmap_d[name] = cmap
    mcm.cmap_d[name + '_r'] = cmap.reversed()
    if re.search('[A-Z]',name):
        mcm.cmap_d[name.lower()] = cmap
        mcm.cmap_d[name.lower() + '_r'] = cmap.reversed()

    # Optionally save colormap to disk
    if name and save:
        if not os.path.isdir(_data_user):
            os.mkdir(_data_user)
        # Save listed colormap i.e. color cycle
        if isinstance(cmap, mcolors.ListedColormap):
            basename = f'{cmap.name}.hex'
            filename = os.path.join(_data_user, basename)
            with open(filename, 'w') as f:
                f.write(','.join(mcolors.to_hex(color) for color in cmap.colors))
        # Save segment data directly
        else:
            basename = f'{cmap.name}.json'
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

def colors(*args, **kwargs):
    """
    Alias for `Cycle`.
    """
    return Cycle(*args, **kwargs)

def Cycle(*args, samples=10, vmin=0, vmax=1, **kwargs):
    """
    Convenience function to draw colors from arbitrary color cycles and
    colormaps stored in `matplotlib.cm.cmap_d`.

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
    # to get samples from a LinearSegmentedColormap
    # draw colors.
    if isinstance(args[-1], Number) or (np.iterable(args[-1])
            and not isinstance(args[-1], (str, dict))):
        args, samples = args[:-1], args[-1]
    # 2) User inputs a simple list; 99% of time, use this
    # to build up a simple ListedColormap.
    elif len(args)>1:
        args = [args] # presumably send a list of colors
    cmap = Colormap(*args, **kwargs) # the cmap object itself
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
    return colors

class PerceptuallyUniformColormap(mcolors.LinearSegmentedColormap):
    """
    Generate a `~matplotlib.colors.LinearSegmentedColormap` in a
    *perceptually uniform colorspace* -- that is, either HSLuv, HCL, or HPLuv.
    See `this page <http://www.hsluv.org/comparison/>`_ for a description
    of each colorspace.
    """
    def __init__(self, name, segmentdata,
            space='hsl', gamma=None, gamma1=None, gamma2=None,
            mask=False, **kwargs):
        """
        Parameters
        ----------
        segmentdata : dict-like
            Dictionary mapping containing the keys ``'hue'``, ``'saturation'``,
            and ``'luminance'``. Values should be lists containing any of
            the following channel specifiers:

                1. Numbers, within the range 0-360 for hue and 0-100 for
                   saturation and luminance.
                2. Color string names or hex tags, in which case the channel
                   value for that color is looked up. See Example.

            See `~matplotlib.colors.LinearSegmentedColormap` for details.
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

        Example
        -------
        The following is a valid `segmentdata` dictionary, using color string
        names for the hue instead of numbers between 0 and 100.

        .. code-block:: python

            dict(hue       = [[0, 'red', 'red'], [1, 'blue', 'blue']],
                saturation = [[0, 1, 1], [1, 1, 1]],
                luminance  = [[0, 1, 1], [1, 0.2, 0.2]])

        Notes
        -----
        Why does `gamma1` emphasize *low* saturation colors, but `gamma2`
        emphasizes *high* luminance colors? Because this seems to be what
        the ColorBrewer2.0 maps do.

        "White" and "pale" values at the
        center of diverging maps and on the left of sequential maps are given
        extra emphasis; the transitions don't seem to be linear in any
        HSL space.
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
        segmentdata['gamma1'] = _default(gamma1, _default(segmentdata.get('gamma1', None), 1.0))
        segmentdata['gamma2'] = _default(gamma2, _default(segmentdata.get('gamma2', None), 1.0))
        self.space = space
        self.mask  = mask
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
                    xyy[j+1] = get_channel_value(y, key, space)
                segmentdata[key][i] = xyy
        # Initialize
        # NOTE: Our gamma1 and gamma2 scaling is just fancy per-channel
        # gamma scaling, so disable the standard version.
        super().__init__(name, segmentdata, gamma=1.0, **kwargs)

    def reversed(self, name=None):
        """
        Returns reversed colormap.
        """
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
        return PerceptuallyUniformColormap(name, data_r, space=self.space)

    def _init(self):
        """
        As with `~matplotlib.colors.LinearSegmentedColormap`, but converts
        each value in the lookup table from 'input' to RGB.
        """
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
            self._lut[i,:3] = to_rgb(self._lut[i,:3], self.space)
        self._lut[:,:3] = clip_colors(self._lut[:,:3], self.mask)

    def _resample(self, N):
        """
        Returns a new colormap with *N* entries.
        """
        self.N = N # that easy
        self._i_under = self.N
        self._i_over = self.N + 1
        self._i_bad = self.N + 2
        self._init()
        return self

    @staticmethod
    def from_hsl(name, h=0, s=100, l=[100, 20], c=None, a=None,
            hue=None, saturation=None, luminance=None, chroma=None, alpha=None,
            ratios=None, reverse=False, **kwargs):
        """
        Make linear segmented colormap by specifying channel values.
        Usage is similar to `~PerceptuallyUniformColormap.from_list` -- except
        instead of specifying HSL tuples, we specify the hue, saturation, and
        luminance vectors or scalars individually.

        Parameters
        ----------
        h : float, str, or list thereof, optional
            Hue channel value or list of values. Values can be
            any of the following:

                1. Numbers, within the range 0-360 for hue and 0-100 for
                   saturation and luminance.
                2. Color string names or hex tags, in which case the channel
                   value for that color is looked up.

            If scalar, the hue does not change across the colormap.
        s : float, str, or list thereof, optional
            As with hue, but for the saturation channel.
        l : float, str, or list thereof, optional
            As with hue, but for the luminance channel.
        a : float, str, or list thereof, optional
            As with hue, but for the alpha channel (the transparency).
        ratios : None or list of float, optional
            Relative extent of each transitions indicated by the channel
            value lists.

            For example, ``luminance=[100,50,0]`` with ``ratios=[2,1]``
            places the *x*-coordinate where the luminance is 50 at 0.66 --
            the white to gray transition is slower than the gray to black
            transition.
        reverse : bool, optional
            Whether to reverse the final colormap.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.

        Other parameters
        ----------------
        hue
            Alias for `h`.
        c, chroma, saturation
            Aliases for `s`.
        luminance
            Alias for `l`.
        alpha
            Alias for `a`.
        """
        # Build dictionary, easy peasy
        h = _default(hue, h)
        s = _default(chroma, _default(c, _default(saturation, s)))
        l = _default(luminance, l)
        a = _default(alpha, _default(a, 1.0))
        cs = ['hue', 'saturation', 'luminance', 'alpha']
        channels = [h, s, l, a]
        cdict = {}
        for c,channel in zip(cs,channels):
            cdict[c] = make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

    @staticmethod
    def from_list(name, color_list,
            ratios=None, reverse=False,
            **kwargs):
        """
        Make linear segmented colormap from list of color tuples.

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
            transition, smaller a faster transition.
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
            cdict[c] = make_segmentdata_array(channel, ratios, reverse, **kwargs)
        cmap = PerceptuallyUniformColormap(name, cdict, **kwargs)
        return cmap

def make_segmentdata_array(values, ratios=None, reverse=False, **kwargs):
    """
    Construct a list of linear segments for an individual channel.
    This was made so that user can input e.g. a callable function for
    one channel, but request linear interpolation for another one.
    """
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
# Colormap constructors
#------------------------------------------------------------------------------#
def merge_cmaps(*imaps, name='merged', N=512, ratios=1, **kwargs):
    """
    Merges arbitrary colormaps. This is used when you pass multiple `args`
    to the `Colormap` function.

    Parameters
    ----------
    *imaps : str or `~matplotlib.colors.Colormap` instances.
        The colormaps for merging.
    name : str, optional
        Name of output colormap. Default is ``'merged'``.
    N : float
        Number of lookup table colors desired for output colormap.

    Notes
    -----
    * The old method had us simply calling the colormap with arrays of fractions.
      This was sloppy, because it just samples locations on the lookup table and
      will therefore degrade the original, smooth, functional transitions.
    * Better method is to combine the ``_segmentdata`` arrays and simply scale
      the *x* coordinates in each ``(x,y1,y2)`` channel-tuple according to
      the ratios.
    * In the case of `~matplotlib.colors.ListedColormap`, we just combine the
      colors.
    """
    # Initial
    if len(imaps)<=1:
        raise ValueError('Need two or more input cmaps.')
    ratios = ratios or 1
    if isinstance(ratios, Number):
        ratios = [1]*len(imaps)

    # Combine the colors
    imaps = [Colormap(cmap, N=None, **kwargs) for cmap in imaps] # set N=None to disable resamping
    if all(isinstance(cmap, mcolors.ListedColormap) for cmap in imaps):
        if not np.all(ratios==1):
            raise ValueError(f'Cannot assign different ratios when mering ListedColormaps.')
        colors = [color for cmap in imaps for color in cmap.colors]
        cmap = mcolors.ListedColormap(colors, name=name, N=len(colors))

    # Accurate methods for cmaps with continuous/functional transitions
    elif all(isinstance(cmap,mcolors.LinearSegmentedColormap) for cmap in imaps):
        # Combine the actual segmentdata
        kinds = {type(cmap) for cmap in imaps}
        if len(kinds)>1:
            raise ValueError(f'Got mixed colormap types.')
        kind = kinds.pop() # colormap kind
        keys = {key for cmap in imaps for key in cmap._segmentdata.keys()}
        ratios = np.array(ratios)/np.sum(ratios) # so if 4 cmaps, will be 1/4
        x0 = np.concatenate([[0], np.cumsum(ratios)])
        xw = x0[1:] - x0[:-1]

        # Combine the segmentdata, and use the y1/y2 slots at merge points
        # so the transition is immediate (can never interpolate between end
        # colors on the two colormaps)
        segmentdata = {}
        for key in keys:
            # Combine scalar values
            if key in ('gamma1', 'gamma2'):
                if key not in segmentdata:
                    segmentdata[key] = []
                for cmap in imaps:
                    segmentdata[key] += [cmap._segmentdata[key]]
                continue
            # Combine xyy data
            datas = []
            test = [callable(cmap._segmentdata[key]) for cmap in imaps]
            if not all(test) and any(test):
                raise ValueError('Mixed callable and non-callable colormap values.')
            if all(test): # expand range from x-to-w to 0-1
                for x,w,cmap in zip(x0[:-1], xw, imaps):
                    data = lambda x: data((x - x0)/w) # WARNING: untested!
                    datas.append(data)
                def data(x):
                    idx, = np.where(x<x0)
                    if idx.size==0:
                        i = 0
                    elif idx.size==x0.size:
                        i = x0.size-2
                    else:
                        i = idx[-1]
                    return datas[i](x)
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

        # Create object
        kwargs = {}
        if kind is PerceptuallyUniformColormap:
            spaces = {cmap.space for cmap in imaps}
            if len(spaces)>1:
                raise ValueError(f'Trying to merge colormaps with different HSL spaces {repr(spaces)}.')
            kwargs.update({'space':spaces.pop()})
        cmap = kind(name, segmentdata, N=N, **kwargs)
    else:
        raise ValueError('All colormaps should be of the same type (Listed or LinearSegmented).')
    return cmap

def monochrome_cmap(color, fade, reverse=False, space='hsl', name='monochrome', **kwargs):
    """
    Make a sequential colormap that blends from some color to near-white.
        Build colormap by varying the luminance of some RGB color while
        keeping its saturation and hue constant.

    Parameters
    ----------
    color : color-like
        Color RGB tuple, hex string, or named color string.
    fade : float or color-like
        The luminance channel strength, or color name from which to take the luminance channel.
    reverse : bool
        Whether to reverse the colormap.
    space : ('hsl')
        Colorspace in which the luminance is varied.
    name : str, optional
        Colormap name. Default is ``'monochrome'``.

    Other parameters
    ----------------
    **kwargs
        Passed to `PerceptuallyUniformColormap.from_hsl` static method.

    Todo
    ----
    Since it's a monochrome colormap, doesn't the HSL colorspace not matter?
    """
    # Get colorspace
    space = _get_space(space)
    h, s, l = to_xyz(to_rgb(color), space)
    if isinstance(fade, Number): # allow just specifying the luminance channel
        fs, fl = s, fade # fade to *same* saturation by default
        fs = s/2
    else:
        _, fs, fl = to_xyz(to_rgb(fade), space)
    index = slice(None,None,-1) if reverse else slice(None)
    return PerceptuallyUniformColormap.from_hsl(name, h,
            [fs,s][index], [fl,l][index], space=space, **kwargs)

def clip_colors(colors, mask=True, gray=0.2, verbose=False):
    """
    Clip impossible colors rendered in an HSl-to-RGB colorspace conversion.

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

    Notes
    -----
    Could use `numpy.clip` (`matplotlib.colors` uses this under the hood),
    but we want to display messages. And anyway premature efficiency is
    the root of all evil, we're manipulating like 1000 colors max here,
    it's no big deal.
    """
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
    # return colors.tolist() # so it is *hashable*, can be cached (wrote this because had weird error, was unrelated)

#------------------------------------------------------------------------------#
# Cycle helper functions
#------------------------------------------------------------------------------#
def set_cycle(cmap, samples=None, rename=False):
    """
    Set the default color cycler.

    Parameters
    ----------
    cmap : str or `~matplotlib.colors.Colormap`
        Colormap from which we draw list of colors.
    samples : None or array-like, optional
        Array of values from 0-1 or number indicating number of evenly spaced
        samples from 0-1 from which to draw colormap colors.

        Will be ignored if the colormap is a `~matplotlib.colors.ListedColormap`,
        for which interpolation is not possible.
    """
    _colors = Cycle(cmap, samples)
    cyl = cycler.cycler('color', _colors)
    rcParams['axes.prop_cycle'] = cyl
    rcParams['patch.facecolor'] = _colors[0]
    if rename:
        rename_colors(cmap)

def rename_colors(cycle='colorblind'):
    """
    Calling this will change how shorthand codes like "b" or "g"
    are interpreted by matplotlib in subsequent plots.

    Parameters
    ----------
    cycle : {'colorblind', 'deep', 'muted', 'bright'}
        The named seaborn palette to use as the source of colors.
    """
    seaborn_cycles = ['colorblind', 'deep', 'muted', 'bright']
    if cycle=='reset':
        colors = [(0.0, 0.0, 1.0), (0.0, .50, 0.0), (1.0, 0.0, 0.0), (.75, .75, 0.0),
                  (.75, .75, 0.0), (0.0, .75, .75), (0.0, 0.0, 0.0)]
    elif cycle in seaborn_cycles:
        colors = cycles[cycle] + [(0.1, 0.1, 0.1)]
    else:
        raise ValueError(f'Cannot set colors with color cycle {cycle}.')
    for code, color in zip('bgrmyck', colors):
        rgb = mcolors.colorConverter.to_rgb(color)
        mcolors.colorConverter.colors[code] = rgb
        mcolors.colorConverter.cache[code]  = rgb

#------------------------------------------------------------------------------#
# Return arbitrary normalizer
#------------------------------------------------------------------------------
def Norm(norm_in, levels=None, values=None, norm=None, **kwargs):
    """
    Return arbitrary normalizer.

    Parameters
    ----------
    norm_in : str or `~matplotlib.colors.Normalize`
        Key name for the normalizer.
    levels, values : array-like
        The level edges (`levels`) or centers (`values`) passed
        to `LinearSegmentedNorm`.
    norm : None or normalizer spec, optional
        The normalizer for *pre-processing*. Used only for
        the `LinearSegmentedNorm` normalizer.

    Other parameters
    ----------------
    **kwargs
        Passed to the `~matplotlib.colors.Normalizer` initializer.
        See `this tutorial <https://matplotlib.org/tutorials/colors/colormapnorms.html>`_
        for more info.

    Notes
    -----
    The recognized normalizer key names are as follows:

    ===============  =================================
    Key              Class
    ===============  =================================
    ``'none'``       `~matplotlib.colors.NoNorm`
    ``'null'``       `~matplotlib.colors.NoNorm`
    ``'zero'``       `MidpointNorm`
    ``'midpoint'``   `MidpointNorm`
    ``'segments'``   `LinearSegmentedNorm`
    ``'segmented'``  `LinearSegmentedNorm`
    ``'boundary'``   `~matplotlib.colors.BoundaryNorm`
    ``'log'``        `~matplotlib.colors.LogNorm`
    ``'linear'``     `~matplotlib.colors.Normalize`
    ``'power'``      `~matplotlib.colors.PowerNorm`
    ``'symlog'``     `~matplotlib.colors.SymLogNorm`
    ===============  =================================

    """
    norm, norm_preprocess = norm_in, norm
    if isinstance(norm, mcolors.Normalize):
        return norm
    if levels is None and values is not None:
        levels = utils.edges(values)
    if not norm: # is None
        # By default, make arbitrary monotonic user levels proceed linearly
        # through color space
        if levels is not None:
            norm = 'segments'
        # Fall back if no levels provided
        else:
            norm = 'linear'
    if isinstance(norm, str):
        # Get class
        if norm not in normalizers:
            raise ValueError(f'Unknown normalizer "{norm}". Options are {", ".join(normalizers.keys())}.')
        norm_out = normalizers[norm]
        # Instantiate class
        if norm_out is BinNorm:
            raise ValueError('This normalizer can only be used internally!')
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
#      generated by your normalizer. This is what BoundaryNorm does.
# Have found the second method was easier to implement/more flexible. So the
# below is used to always discretize colors.
#------------------------------------------------------------------------------
# WARNING: Many methods in ColorBarBase tests for class membership, crucially
# including _process_values(), which if it doesn't detect BoundaryNorm will
# end up trying to infer boundaries from inverse() method. So make it parent class.
class BinNorm(mcolors.BoundaryNorm):
    """
    This is a rough copy of BoundaryNorm, but includes some extra features.
    *Discreteizes* the possible normalized values (numbers in 0-1) that are
    used to index a color on a high-resolution colormap lookup table. But
    includes features for other stuff.

    Note
    ----
    If you are using a diverging colormap with extend='max/min', the center
    will get messed up. But that is very strange usage anyway... so please
    just don't do that :)

    Todo
    ----
    Allow this to accept transforms too, which will help prevent level edges
    from being skewed toward left or right in case of logarithmic/exponential data.

    Example
    -------
    Your levels edges are weirdly spaced [-1000, 100, 0, 100, 1000] or
    even [0, 10, 12, 20, 22], but center "colors" are always at colormap
    coordinates [.2, .4, .6, .8] no matter the spacing; levels just must be monotonic.
    """
    def __init__(self, levels, norm=None, clip=False, step=1.0, extend='neither', **kwargs):
        # Declare boundaries, vmin, vmax in True coordinates. The step controls
        # intensity transition to out-of-bounds color; by default, the step is
        # equal to the *average* step between in-bounds colors (step == 1).
        # NOTE: Idea is that we bin data into len(levels) discrete x-coordinates,
        # and optionally make out-of-bounds colors the same or different
        # NOTE: Don't need to call parent __init__, this is own implementation
        # Do need it to subclass BoundaryNorm, so ColorbarBase will detect it
        # See BoundaryNorm: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        extend = extend or 'both'
        levels = np.atleast_1d(levels)
        if levels.size<=1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError(f'Levels {levels} passed to Normalize() must be monotonically increasing.')
        if extend not in ('both','min','max','neither'):
            raise ValueError(f'Unknown extend option "{extend}". Choose from "min", "max", "both", "neither".')

        # Determine color ids for levels, i.e. position in 0-1 space
        # NOTE: If user used LinearSegmentedNorm for the normalizer (the
        # default) any monotonic levels will be even.
        # NOTE: Length of these ids should be N + 1 -- that is, N - 1 colors
        # for values in-between levels, plus 2 colors for out-of-bounds.
        #   * For same out-of-bounds colors, looks like [0, 0, ..., 1, 1]
        #   * For unique out-of-bounds colors, looks like [0, X, ..., 1 - X, 1]
        #     where the offset X equals step/len(levels).
        # First get coordinates
        norm = norm or (lambda x: x) # e.g. a logarithmic transform
        x_b = norm(levels)
        x_m = (x_b[1:] + x_b[:-1])/2 # get level centers after norm scaling
        y = (x_m - x_m.min())/(x_m.max() - x_m.min())
        # Account for out of bounds colors
        offset = 0
        scale = 1
        eps = step/levels.size
        if extend in ('min','both'):
            offset = eps
            scale -= eps
        if extend in ('max','both'):
            scale -= eps
        if isinstance(y, ma.core.MaskedArray):
            y = y.filled(np.nan)
        y = y[np.isfinite(y)]
        y = np.concatenate(([0], offset + scale*y, [1])) # insert '0' (arg 3) before index '0' (arg 2)
        self._norm = norm
        self._x_b = x_b
        self._y = y

        # Add builtin properties
        # NOTE: Are vmin/vmax even used?
        self.boundaries = levels
        self.vmin = levels.min()
        self.vmax = levels.max()
        self.clip = clip
        self.N = levels.size

    def __call__(self, xq, clip=None):
        """
        Normalize data values to the range 0-1.
        """
        # Follow example of LinearSegmentedNorm, but perform no interpolation,
        # just use searchsorted to bin the data.
        # NOTE: The bins vector includes out-of-bounds negative (searchsorted
        # index 0) and out-of-bounds positive (searchsorted index N+1) values
        xq = self._norm(np.atleast_1d(xq))
        yq = self._y[np.searchsorted(self._x_b, xq)] # which x-bin does each point in xq belong to?
        return ma.masked_array(yq, np.isnan(xq))

    def inverse(self, yq):
        """
        Dummy method.
        """
        # Not possible
        raise ValueError('BinNorm is not invertible.')

#------------------------------------------------------------------------------#
# Normalizers intended to *pre-scale* levels passed to BinNorm
#------------------------------------------------------------------------------#
class LinearSegmentedNorm(mcolors.Normalize):
    """
    Linearly *interpolate* colors between the provided boundary levels.
    Exactly analagous to the method in `~matplotlib.colors.LinearSegmentedColormap`:
    performs linear interpolation between successive monotonic, but arbitrarily
    spaced, points. Then lets this index control color intensity.

    This class is useful when you want "evenly spaced" colors for unevenly
    spaced levels -- e.g. when your data spans a large range of magnitudes.

    This class is the **default** normalizer for all functions that accept
    the `cmap` keyword arg.
    """
    def __init__(self, levels, clip=False, **kwargs):
        # Test
        levels = np.atleast_1d(levels)
        if levels.size<=1:
            raise ValueError('Need at least two levels.')
        elif ((levels[1:]-levels[:-1])<=0).any():
            raise ValueError(f'Levels {levels} passed to LinearSegmentedNorm must be monotonically increasing.')
        super().__init__(np.nanmin(levels), np.nanmax(levels), clip) # second level superclass
        self._x = levels
        self._y = np.linspace(0, 1, len(levels))

    def __call__(self, xq, clip=None):
        """
        Normalize data values to the range 0-1.
        """
        # Follow example of make_mapping_array for efficient, vectorized
        # linear interpolation across multiple segments
        # NOTE: normal test puts values at a[i] if a[i-1] < v <= a[i]; for
        # left-most data, satisfy a[0] <= v <= a[1]
        # NOTE: searchsorted gives where xq[i] must be inserted so it is larger
        # than x[ind[i]-1] but smaller than x[ind[i]]
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
        """
        Inverse operation of `__call__`.
        """
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
    Ensure a "midpoint" always lies at the central colormap color. This
    is normally used with diverging colormaps and ``midpoint=0``.

    Notes
    -----
    See `this stackoverflow thread <https://stackoverflow.com/q/25500541/4970632>`_.
    """
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=None):
        # Bigger numbers are too one-sided
        super().__init__(vmin, vmax, clip)
        self._midpoint = midpoint

    def __call__(self, xq, clip=None):
        """
        Normalize data values to the range 0-1.
        """
        # Get middle point in 0-1 coords, and value
        # NOTE: Look up these three values in case vmin/vmax changed; this is
        # a more general normalizer than the others. Others are 'parent'
        # normalizers, meant to be static more or less.
        # NOTE: searchsorted gives where xq[i] must be inserted so it is larger
        # than x[ind[i]-1] but smaller than x[ind[i]]
        # x, y = [self.vmin, self._midpoint, self.vmax], [0, 0.5, 1]
        # return ma.masked_array(np.interp(xq, x, y))
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

    def inverse(self, yq, clip=None):
        """
        Inverse operation of `__call__`.
        """
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
    Register new color names and **filter** them to be necessarily
    "perceptually distinct" in the HSL colorspace.

    Use `~proplot.demos.color_show` to generate a table of the resulting
    filtered colors.

    Notes
    -----
    This seems like it would be slow, but takes on average 0.03 seconds on
    my macbook, so it's fine.
    """
    # First ***reset*** the colors dictionary
    # Why? We want to add XKCD colors *sorted by popularity* from file, along
    # with crayons dictionary; having registered colors named 'xkcd:color' is
    # annoying and not useful
    scale = (360, 100, 100)
    translate =  {'b': 'blue', 'g': 'green', 'r': 'red', 'c': 'cyan',
                  'm': 'magenta', 'y': 'yellow', 'k': 'black', 'w': 'white'}
    base1 = mcolors.BASE_COLORS # one-character names
    base2 = {translate[key]:value for key,value in base1.items()} # full names
    mcolors._colors_full_map.clear() # clean out!
    mcolors._colors_full_map.cache.clear() # clean out!
    mcolors._colors_full_map.update(base1)
    mcolors._colors_full_map.update(base2)

    # First register colors and get their HSL values
    # Sort files in reverse order becuase I prefer XKCD color names to crayon
    # color names; so want to overwrite identical names with XKCD name.
    seen = set()
    names = []
    hcls = np.empty((0,3))
    correct = (('/', ' '), ("'s", ''), ('grey', 'gray'),
               ('pinky', 'pink'), ('greeny', 'green'),
               ('robin egg', 'robins egg'),
               ('egg blue', 'egg'),
               (r'reddish', 'red'),
               (r'purplish', 'purple'),
               (r'bluish',  'blue'),
               (r'ish\b', ''),
               ('bluegray', 'blue gray'),
               ('grayblue', 'gray blue'),
               ('lightblue', 'light blue'))
    for file in sorted(glob.glob(os.path.join(_data_colors, '*.txt'))):
        # Read data
        category, _ = os.path.splitext(os.path.basename(file))
        data = np.genfromtxt(file, delimiter='\t', dtype=str, comments='%', usecols=(0,1)).tolist()
        # Sanitize names and add to dictionary
        # NOTE: You can add to this!
        i = 0
        _dict = {}
        ihcls = []
        colorlist[category] = {} # just initialize this one
        for name, color in data: # is list of name, color tuples
            if i>=nmax: # e.g. for xkcd colors
                break
            # Sanitize
            for reg,sub in correct:
                name = re.sub(reg, sub, name)
            # Check
            if name in seen:
                continue
            # Add
            seen.add(name)
            names.append((category, name)) # save the category name pair
            ihcls.append(to_xyz(color, space=_distinct_colors_space))
            _dict[name] = color # save the color
            i += 1
        _colors_unfiltered[category] = _dict
        # Concatenate HCL arrays
        hcls = np.concatenate((hcls, ihcls), axis=0)

    # Remove colors that are 'too similar' by rounding to the nearest n units
    # WARNING: unique axis argument requires numpy version >=1.13
    # WARNING: evidently it is ***impossible*** to actually delete colors
    # from the custom_colors dictionary (perhaps due to quirk of autoreload,
    # perhaps by some more fundamental python thing), so we instead must create
    # *completely separate* dictionary and add colors from there
    hcls = hcls/np.array(scale)
    hcls = np.round(hcls/_distinct_colors_threshold).astype(np.int64)
    _, index, counts = np.unique(hcls, return_index=True, return_counts=True, axis=0) # get unique rows
    deleted = 0
    counts = counts.sum()
    exceptions_regex = '^(' + '|'.join(_distinct_colors_exceptions) + ')[0-9]?$'

    # Add colors to filtered colors
    for i,(category,name) in enumerate(names):
        if not re.match(exceptions_regex, name) and i not in index:
            deleted += 1
        else:
            colorlist[category][name] = _colors_unfiltered[category][name]
    for key,kw in colorlist.items():
        mcolors._colors_full_map.update(kw)
    if verbose:
        print(f'Started with {len(names)} colors, removed {deleted} insufficiently distinct colors.')

def register_cmaps():
    """
    Register colormaps and color cycles in the `cmaps` directory packaged
    with ProPlot. That is, add colors to the `matplotlib.cm.cmap_d`
    dictionary.

    Use `~proplot.demos.cmap_show` to generate a table of the resulting
    color cycles.
    """
    # First read from file
    for filename in sorted(glob.glob(os.path.join(_data_cmaps, '*'))) + \
            sorted(glob.glob(os.path.join(_data_user, '*'))):
        # Read table of RGB values
        if not re.search('\.(x?rgba?|json|xml)$', filename):
            continue
        name = os.path.basename(filename)
        name = name.split('.')[0]
        # if name in mcm.cmap_d: # don't want to re-register every time
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
                cmap = PerceptuallyUniformColormap(name, segmentdata, space=space, N=_N_hires)
            else:
                cmap = mcolors.LinearSegmentedColormap(name, segmentdata, N=_N_hires)
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
    for name in ['twilight', 'Phase']:
        cmap = mcm.cmap_d.get(name, None)
        if cmap and isinstance(cmap, mcolors.LinearSegmentedColormap):
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
            mcm.cmap_d[name] = mcolors.LinearSegmentedColormap(name, data, cmap.N)
            mcm.cmap_d[name + '_shifted'] = mcolors.LinearSegmentedColormap(name + '_shifted', data_shift, cmap.N)

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
    Register color cycles defined in ``.hex`` files (lists of hex strings)
    and those declared in this module.

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
        colors = [*open(filename)] # should just be a single line
        if len(colors)==0:
            continue # file is empty
        if len(colors)>1:
            raise ValueError('.hex color cycle files should contain only one line.')
        colors = colors[0].strip().split(',') # csv hex strings
        colors = [mcolors.to_rgb(c) for c in colors] # from list of tuples
        _cycles_loaded[name] = colors

    # Register names
    # Note that 'loaded' cycles will overwrite any presets with same name
    for name,colors in {**_cycles_preset, **_cycles_loaded}.items():
        mcm.cmap_d[name] = mcolors.ListedColormap([to_rgb(color) for color in colors])
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
colorlist = {} # limit to 'sufficiently unique' color names
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
    'boundary':   mcolors.BoundaryNorm,
    'log':        mcolors.LogNorm,
    'linear':     mcolors.Normalize,
    'power':      mcolors.PowerNorm,
    'symlog':     mcolors.SymLogNorm,
    }
"""Dictionary of possible normalizers."""

