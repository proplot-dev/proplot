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

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

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
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # 'sphinx.ext.ifconfig',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon', # for NumPy style docstrings, instead of reStructred Text
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autosummary',
    'sphinxcontrib.bibtex', # see: https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html
    'sphinxext.automodapi', # see: https://sphinx-automodapi.readthedocs.io/en/latest/
    # 'nbsphinx',
    # 'IPython.sphinxext.ipython_directive', # for ipython highlighting
    # 'IPython.sphinxext.ipython_console_highlighting',
    # 'sphinxext.custom_roles', # copied directly from matplotlib
    # 'matplotlib.sphinxext.only_directives', # deprecated; see: https://github.com/statsmodels/statsmodels/issues/5291
    # 'matplotlib.sphinxext.plot_directive', # see: https://matplotlib.org/sampledoc/extensions.html
    ]

# Generate stub pages whenever ::autosummary directive encountered
# This way don't have to call sphinx-autogen manually
autosummary_generate = True

# Use automodapi tool, created by astropy people
# See: https://sphinx-automodapi.readthedocs.io/en/latest/automodapi.html#overview
# Normally have to *enumerate* function names manually. This will document
# them automatically. Just be careful, if you use from x import *, to exclude
# them in the automodapi:: directive
# Also modify so 
automodapi_toctreedirnm = 'api' # create much better URL for the page
automodsumm_inherited_members = False
automodsumm_inherited_members = False

# Logo
html_logo = '_static/logo.png'

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
napoleon_use_rtype = False
napoleon_use_param = False
napoleon_use_ivar = False
napoleon_use_keyword = False
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False # move init doc to 'class' doc

# Fix duplicate class member documentation from autosummary + numpydoc
# See: https://github.com/phn/pytpm/issues/3#issuecomment-12133978
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
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
exclude_patterns = ['_build', '_templates', '_themes', 'sphinxext',
                    'originals', 'showcase', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# Role
default_role = 'py:obj' # default family is py, but can also set default role so don't need :func:`name`, :module:`name`, etc.

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster' # ugly default
# html_theme = 'sphinx_rtd_theme' # from: https://sphinx-rtd-theme.readthedocs.io/en/latest/installing.html
html_theme = 'rtd_custom'
html_theme_path = ['_themes']
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    }
# html_sidebars = {'**': ['layout.html', 'search.html', 'breadcrumbs.html', 'footer.html', 'searchbox.html']}
# html_sidebars = {'**': ['fulltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']}
# html_sidebars = {'**': ['fulltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
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
# To convert: convert icon.png icon.ico
html_favicon = os.path.join('_static', 'icon.ico')

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
