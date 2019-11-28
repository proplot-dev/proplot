# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# For autodoc compilation see:
# https://medium.com/@eikonomega/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

import os
import sys
import matplotlib # load local matplotlibrc and set up docstring settings
from pygments.formatters import HtmlFormatter
from pygments.styles import get_all_styles

# Sphinx-automodapi requires proplot on path
sys.path.insert(0, os.path.abspath('..'))
# Then add path for local 'sphinxext' extensions
# Not sure when abspath is required
sys.path.append(os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'ProPlot'
copyright = '2019, Luke L. B. Davis'
author = 'Luke L. B. Davis'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    # 'matplotlib.sphinxext.plot_directive',  # see: https://matplotlib.org/sampledoc/extensions.html
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',  # for ipython highlighting
    'sphinx.ext.autodoc',           # include documentation from docstrings
    'sphinx.ext.doctest',           # >>> examples
    'sphinx.ext.extlinks',          # for :pr:, :issue:, :commit:
    'sphinx.ext.autosectionlabel',  # use :ref:`Heading` for any heading
    'sphinx.ext.todo',              # Todo headers and todo:: directives
    'sphinx.ext.mathjax',           # LaTeX style math
    'sphinx.ext.viewcode',          # view code links
    'sphinx.ext.autosummary',       # autosummary directive
    'sphinx.ext.napoleon',          # for NumPy style docstrings, instead of reStructred Text
    'sphinx.ext.intersphinx',       # external links
    'sphinxext.custom_roles',       # local extension
    'sphinx_automodapi.automodapi', # see: https://github.com/lukelbd/sphinx-automodapi/tree/proplot-mods
    'nbsphinx',
    ]

extlinks = {
    'issue': ('https://github.com/lukelbd/proplot/issues/%s', 'GH#'),
    'commit': ('https://github.com/lukelbd/proplot/commit/%s', '@'),
    'pr': ('https://github.com/lukelbd/proplot/pull/%s', 'GH#'),
}

# Give *lots* of time for cell execution! The projection tables
# in particular are massive.
nbsphinx_timeout = 120

# Do not run doctest tests, these are just to show syntax and expected
# output may be graphical
doctest_test_doctest_blocks = ''

# Generate stub pages whenever ::autosummary directive encountered
# This way don't have to call sphinx-autogen manually
autosummary_generate = True

# Use automodapi tool, created by astropy people
# See: https://sphinx-automodapi.readthedocs.io/en/latest/automodapi.html#overview
# Normally have to *enumerate* function names manually. This will document
# them automatically. Just be careful, if you use from x import *, to exclude
# them in the automodapi:: directive
automodapi_toctreedirnm = 'api' # create much better URL for the page
automodsumm_inherited_members = False

# Logo
html_logo = '_static/logo_square.png'

# Turn off code and image links for embedded mpl plots
# plot_html_show_source_link = False
# plot_html_show_formats = False

# One of 'class', 'both', or 'init'
# The 'both' concatenates class and __init__ docstring
# See: http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
autoclass_content = 'both'

# Set up mapping for other projects' docs
intersphinx_mapping = {
                       'cycler': ('https://matplotlib.org/cycler/', None),
                       'matplotlib': ('https://matplotlib.org', None),
                       'sphinx': ('http://www.sphinx-doc.org/en/stable', None),
                       'python': ('https://docs.python.org/3', None),
                       'numpy': ('https://docs.scipy.org/doc/numpy', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
                       'xarray': ('http://xarray.pydata.org/en/stable', None),
                       'cartopy': ('https://scitools.org.uk/cartopy/docs/latest', None),
                       'basemap': ('https://matplotlib.org/basemap', None),
                       'pandas': ('https://pandas.pydata.org/pandas-docs/stable', None),
                       }

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False # confusing, because I use submodules for *organization*

# Napoleon options
# See: http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
napoleon_use_ivar = False
napoleon_use_param = False
napoleon_use_keyword = False
napoleon_use_rtype = False
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False # move init doc to 'class' doc

# Fix duplicate class member documentation from autosummary + numpydoc
# See: https://github.com/phn/pytpm/issues/3#issuecomment-12133978
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
# WARNING: Must add 'include' files or will get duplicate label warnings.
# WARNING: Must add files containing showcase examples
exclude_patterns = [
    '_templates', '_themes', 'showcase',
    'sphinxext', 'automodapi',
    'trash', '.DS_Store', '**.ipynb_checkpoints'
    ]

# The name of the Pygments (syntax highlighting) style to use.
# The light-dark theme toggler overloads this, but set default anyway
pygments_style = 'none'

# Create local pygments copies
# Previously used: https://github.com/richleland/pygments-css
# But do not want to depend on some random repository
path = os.path.join('_static', 'pygments')
if not os.path.isdir(path):
    os.mkdir(path)
for style in get_all_styles():
    path = os.path.join('_static', 'pygments', style + '.css')
    if os.path.isfile(path):
        continue
    with open(path, 'w') as f:
        f.write(HtmlFormatter(style=style).get_style_defs('.highlight'))

# Role
default_role = 'py:obj' # default family is py, but can also set default role so don't need :func:`name`, :module:`name`, etc.

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# Rtd theme still the best
# in _templates but can just use below optoin.
# We set "style_nav_header_background" in custom.css
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    'prev_next_buttons_location': 'bottom', # top and bottom
    }

# Matplotlib theme, ugly
# html_theme = "sphinxdoc" # like matplotlib

# Beautiful but sidebar is not locked
# import sphinx_modern_theme
# html_theme = 'sphinx_modern_theme'
# pygments_style = 'default'
# html_theme_path = [
#     sphinx_modern_theme.get_html_theme_path(),
# ]

# Simple but insufficient
# Tried: https://stackoverflow.com/a/57040610/4970632
# But items disappear after adding html_sidebars dictionary
# html_theme = 'alabaster'
# html_logo = None # or get 2 logos! need to specify in dictionary so text does not appear
# html_theme_options = {
#     'logo': 'logo_square.png',
#     'logo_name': False,
#     'description': 'A matplotlib wrapper for making beautiful, publication-quality graphics.',
#     'page_width': '90%',
#     'fixed_sidebar': True,
#     'show_relbars': True,
#     'sidebar_collapse': True,
#     'sidebar_width': '20%',
#     'caption_font_size': 'x-large',
#     'sidebar_link_underscore': 'transparent',
#     }

# Clean but no margins, no full TOC, no box around tables
# To modify see: For guzzle theme: https://github.com/guzzle/guzzle_sphinx_theme/issues/22
# import guzzle_sphinx_theme
# html_theme_path = guzzle_sphinx_theme.html_theme_path()
# html_theme = 'guzzle_sphinx_theme'
# extensions.append("guzzle_sphinx_theme")
# html_theme_options = {
#     # Set the name of the project to appear in the sidebar
#     # "project_nav_name": "Project Name",
#     'navigation_depth': 4,
# }

# Readthedocs clones
# import rtcat_sphinx_theme
# html_theme = "rtcat_sphinx_theme"
# html_theme_path = [rtcat_sphinx_theme.get_html_theme_path()]
# import sphinx_pdj_theme
# html_theme = 'sphinx_pdj_theme'
# htm_theme_path = [sphinx_pdj_theme.get_html_theme_path()]

# Custom theme in the future?
# html_theme = 'custom'
# html_theme_path = ['_themes']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large. Static folder is for CSS and image files.
# For icons see: https://icons8.com/icon
# To convert: convert logo_blank.png logo_blank.ico
html_favicon = os.path.join('_static', 'logo_blank.ico')

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'proplotdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'proplot.tex', 'proplot Documentation',
     'Luke L. B. Davis', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'proplot', 'proplot Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'proplot', 'proplot Documentation',
     author, 'proplot', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True
