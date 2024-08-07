# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts:
#
#
# Font selection
# ==============
#
# Proplot registers several new fonts and includes tools
# for adding your own fonts. These features are described below.
#
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts_included:
#
# Included fonts
# --------------
#
# Matplotlib provides a `~matplotlib.font_manager` module for working with
# system fonts and classifies fonts into `five font families
# <https://matplotlib.org/stable/gallery/text_labels_and_annotations/fonts_demo.html>`__:
# :rcraw:`font.serif` :rcraw:`font.sans-serif`, :rcraw:`font.monospace`,
# :rcraw:`font.cursive`, and :rcraw:`font.fantasy`. The default font family
# is sans-serif, because sans-serif fonts are generally more suitable for
# figures than serif fonts, and the default font name belonging to this family
# is `DejaVu Sans <https://dejavu-fonts.github.io>`__, which comes packaged with
# matplotlib.
#
# Matplotlib uses DejaVu Sans in part because it includes glyphs for a very wide
# range of symbols, especially mathematical symbols. However in our opinion,
# DejaVu Sans is not very aesthetically pleasing. To improve the font selection while
# keeping things consistent across different workstations, proplot is packaged
# the open source `TeX Gyre fonts <https://ctan.org/pkg/tex-gyre?lang=en>`__ and a few
# additional open source sans-serif fonts. Proplot also uses the TeX Gyre fonts as the
# first (i.e., default) entries for each of matplotlib's `font family lists
# <https://matplotlib.org/stable/tutorials/text/text_props.html#default-font>`__:
#
# * The `Helvetica <https://en.wikipedia.org/wiki/Helvetica>`__ lookalike
#   :rcraw:`font.sans-serif` = ``'TeX Gyre Heros'``.
# * The `Century <https://en.wikipedia.org/wiki/Century_type_family>`__ lookalike
#   :rcraw:`font.serif` = ``'TeX Gyre Schola'``.
# * The `Chancery <https://en.wikipedia.org/wiki/ITC_Zapf_Chancery>`__ lookalike
#   :rcraw:`font.cursive` = ``'TeX Gyre Chorus'``.
# * The `Avant Garde <https://en.wikipedia.org/wiki/ITC_Avant_Garde>`__ lookalike
#   :rcraw:`font.fantasy` = ``'TeX Gyre Adventor'``.
# * The `Courier <https://en.wikipedia.org/wiki/Courier_(typeface)>`__ lookalike
#   :rcraw:`font.monospace` = ``'TeX Gyre Cursor'``.
#
# After importing proplot, the default matplotlib font will be
# `TeX Gyre Heros <https://ctan.org/pkg/tex-gyre-heros>`__, which
# emulates the more conventional and (in our opinion) aesthetically pleasing
# font `Helvetica <https://en.wikipedia.org/wiki/Helvetica>`__. The default font
# family lists are shown in the :ref:`default proplotrc file <ug_proplotrc>`.
# To compare different fonts, use the `~proplot.demos.show_fonts` command with the
# `family` keyword (default behavior is ``family='sans-serif'``). Tables of the TeX
# Gyre and sans-serif fonts packaged with proplot are shown below.

# %%
import proplot as pplt

fig, axs = pplt.show_fonts(family="sans-serif")

# %%
import proplot as pplt

fig, axs = pplt.show_fonts(family="tex-gyre")

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts_math:
#
# Math text fonts
# ---------------
#
# In matplotlib, math text rendered by TeX can be produced by surrounding
# an expression with ``$dollar signs$``. To help math text jive better with
# the new default :ref:`non-math text font <ug_fonts_included>`, proplot changes
# :rcraw:`mathtext.fontset` to ``'custom'``. This means that math is drawn with
# the italicized version of the non-math font (see the matplotlib `math text
# guide <https://matplotlib.org/stable/tutorials/text/mathtext.html#custom-fonts>`__
# for details). This generally improves the appearance of figures with simple
# math expressions. However, if you need unusual math symbols or complex math
# operators, you may want to change :rcraw:`font.name` to something more suitable
# for math (e.g., the proplot-packaged font ``'Fira Math'`` or the matplotlib-packaged
# font ``'DejaVu Sans'``; see `this page <https://github.com/firamath/firamath>`__ for
# more on Fira Math). Alternatively, you can change the math text font alone by setting
# :rcraw:`mathtext.fontset` back to one of matplotlib's math-specialized font sets
# (e.g., ``'stixsans'`` or ``'dejavusans'``).
#
# A table of math text containing the sans-serif fonts packaged with proplot is shown
# below. The dummy glyph "¤" is shown where a given math character is unavailable
# for a particular font (in practice, the fallback font :rc:`mathtext.fallback` is used
# whenever a math character is unavailable, but `~proplot.demos.show_fonts` disables
# this fallback font in order to highlight the missing characters).
#
# .. note::
#
#    Proplot modifies matplotlib's math text internals so that the ``'custom'``
#    font set can be applied with modifications to the currently active non-math
#    font rather than only a global font family. This works by changing the default
#    values of :rcraw:`mathtext.bf`, :rcraw:`mathtext.it`, :rcraw:`mathtext.rm`,
#    :rcraw:`mathtext.sf` from the global default font family ``'sans'`` to the local
#    font family ``'regular'``, where ``'regular'`` is a dummy name permitted by
#    proplot (see the :ref:`proplotrc file <ug_proplotrc>` for details). This means
#    that if :rcraw:`mathtext.fontset` is ``'custom'`` and the font family is changed
#    for an arbitrary `~matplotlib.text.Text` instance, then any LaTeX-generated math
#    in the text string will also use this font family.

# %%
import proplot as pplt

fig, axs = pplt.show_fonts(family="sans-serif", math=True)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts_user:
#
# Using your own fonts
# --------------------
#
# You can register your own fonts by adding files to the ``fonts`` subfolder
# inside `~proplot.config.Configurator.user_folder` and calling
# `~proplot.config.register_fonts`. This command is called on import. You can
# also manually pass file paths to `~proplot.config.register_fonts`.
# To change the default font, use the `~proplot.config.rc`
# object or modify your ``proplotrc``. See the
# :ref:`configuration section <ug_config>` for details.
#
# Sometimes the font you would like to use *is* installed, but the font file
# is not stored under the matplotlib-compatible ``.ttf``, ``.otf``, or ``.afm``
# formats. For example, several macOS fonts are unavailable because they are
# stored as ``.dfont`` collections. Also, while matplotlib nominally supports
# ``.ttc`` collections, proplot ignores them because figures with ``.ttc`` fonts
# `cannot be saved as PDFs <https://github.com/matplotlib/matplotlib/issues/3135>`__.
# You can get matplotlib to use ``.dfont`` and ``.ttc`` collections by
# expanding them into individual ``.ttf`` files with the
# `DFontSplitter application <https://peter.upfold.org.uk/projects/dfontsplitter>`__,
# then saving the files in-place or in the ``~/.proplot/fonts`` folder.
#
# To find font collections, check the paths listed in ``OSXFontDirectories``,
# ``X11FontDirectories``, ``MSUserFontDirectories``, and ``MSFontDirectories``
# under the `matplotlib.font_manager` module.
