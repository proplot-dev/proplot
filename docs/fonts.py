# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
#
# Font selection
# ==============
#
# ProPlot registers several new font families and includes tools for adding
# your own fonts. These features are described below.
#
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts:
#
# Included fonts
# --------------
#
# Matplotlib provides a `~matplotlib.font_manager` module for working with
# system fonts and classifies fonts into `five font families\
# <https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/fonts_demo.html>`__:
# :rcraw:`font.serif` :rcraw:`font.sans-serif`, :rcraw:`font.monospace`,
# :rcraw:`font.cursive`, and :rcraw:`font.fantasy`. The default font family
# is sans-serif, because sans-serif fonts are generally more suitable for
# figures than serif fonts, and the default font name belonging to this family
# is `DejaVu Sans <https://dejavu-fonts.github.io>`__, which comes packaged with
# matplotlib.
#
# Matplotlib uses DejaVu Sans in part because it includes glyphs
# for a wider range of mathematical symbols. However in your opinion, DejaVu Sans
# not very aesthetically pleasing or particularly readable, and it is
# seldom used outside of the matplotlib ecosystem.
# To improve the font selection and make things consistent across different
# workspaces, ProPlot comes packaged with the open-source
# `TeX Gyre <http://www.gust.org.pl/projects/e-foundry/tex-gyre>`__ font series
# and adds them as the default entries for all of matplotlib's font famlies:
#
# * The `Century <https://en.wikipedia.org/wiki/Century_type_family>`__ lookalike
#   :rcraw:`font.serif` = ``'TeX Gyre Schola'``.
# * The `Helvetica <https://en.wikipedia.org/wiki/Helvetica>`__ lookalike
#   :rcraw:`font.sans-serif` = ``'TeX Gyre Heros'``.
# * The `Courier <https://en.wikipedia.org/wiki/Courier_(typeface)>`__ lookalike
#   :rcraw:`font.monospace` = ``'TeX Gyre Cursor'``.
# * The `Chancery <https://en.wikipedia.org/wiki/ITC_Zapf_Chancery>`__ lookalike
#   :rcraw:`font.cursive` = ``'TeX Gyre Chorus'``.
# * The `Avant Garde <https://en.wikipedia.org/wiki/ITC_Avant_Garde>`__ lookalike
#   :rcraw:`font.fantasy` = ``'TeX Gyre Adventor'``.
#
# Thus after importing ProPlot, even on sparse Linux servers with limited font
# selection, the default font will be a more conventional and aesthetically
# pleasing Helvetica lookalike. The new font priority lists for each font family
# are shown in the :ref:`default proplotrc file <ug_proplotrc>`.
#
# To compare different fonts, use the `~proplot.show.show_fonts` command. By
# default, this displays the *sans-serif* fonts available on your system and
# included with ProPlot. The default table on a sparse Linux server is shown
# below. The "¤" symbol appears where characters for a particular font are
# unavailable. When you are actually making plots, "¤" will be filled with
# the character from a fallback font.
#
# It should be noted that most of the TeX Gyre fonts have limited character
# availability. If your plots contain lots of mathematical symbols,
# you may want to switch the default :rcraw:`font.family` to DejaVu Sans, which
# is packaged with matplotlib, or `Fira Math <https://github.com/firamath/firamath>`__,
# which is also packaged with ProPlot.

# %%
import proplot as plot
fig = plot.show_fonts()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts_user:
#
# Using your own fonts
# --------------------
#
# You can register your own fonts by adding files to the ``~/.proplot/fonts``
# directory and calling `~proplot.config.register_fonts`. This command is
# also called on import. To change the default font, use the
# `~proplot.config.rc` object or modify your ``~/.proplotrc``. See
# the :ref:`configuration section <ug_config>` for details.
#
# Sometimes the font you would like to use *is* installed, but the font file
# is not stored under the matplotlib-compatible ``.ttf``, ``.otf``, or
# ``.afm`` formats. For example, several macOS fonts are unavailable because
# they are stored as ``.dfont`` collections. Also, while matplotlib nominally
# supports ``.ttc`` collections, ProPlot manually removes them because
# figures with ``.ttc`` fonts
# `cannot be saved as PDFs <https://github.com/matplotlib/matplotlib/issues/3135>`__.
# You can get matplotlib to use these fonts by expanding the "collections"
# into individual ``.ttf`` files with the
# `DFontSplitter application <https://peter.upfold.org.uk/projects/dfontsplitter>`__,
# then saving the files in-place or in the ``~/.proplot/fonts`` folder.
#
# To find font files, check the paths listed in ``OSXFontDirectories``,
# ``X11FontDirectories``, ``MSUserFontDirectories``, and ``MSFontDirectories``
# under the `~matplotlib.font_manager` module. Note that if the font in question has
# a "thin" style, implied by file names with the word ``Thin``,
# `a matplotlib bug <https://github.com/matplotlib/matplotlib/issues/8788>`__
# may cause these styles to override the "normal" style.
