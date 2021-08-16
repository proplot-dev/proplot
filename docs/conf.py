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
import datetime

# Add proplot to path for sphinx-automodapi
sys.path.insert(0, os.path.abspath('..'))

# Add docs folder to PATH for local 'sphinxext' extensions
sys.path.append(os.path.abspath('.'))

# Hack to get basemap to work
# See: https://github.com/readthedocs/readthedocs.org/issues/5339
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    os.environ['PROJ_LIB'] = (
        '{}/{}/share/proj'.format(
            os.environ['CONDA_ENVS_PATH'], os.environ['CONDA_DEFAULT_ENV']
        )
    )
else:
    os.environ['PROJ_LIB'] = '{}/share/proj'.format(
        os.environ['CONDA_PREFIX']
    )

# -- Project information -----------------------------------------------------

_today = datetime.datetime.today()
project = 'ProPlot'
copyright = f'{_today.year}, Luke L. B. Davis'
author = 'Luke L. B. Davis'

# The short X.Y version
version = ''

# The full version, including alpha/beta/rc tags
release = ''

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    # 'matplotlib.sphinxext.plot_directive', # see: https://matplotlib.org/sampledoc/extensions.html  # noqa
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',  # for ipython highlighting
    'sphinx.ext.autodoc',  # include documentation from docstrings
    'sphinx.ext.doctest',  # >>> examples
    'sphinx.ext.extlinks',  # for :pr:, :issue:, :commit:
    'sphinx.ext.autosectionlabel',  # use :ref:`Heading` for any heading
    'sphinx.ext.todo',  # Todo headers and todo:: directives
    'sphinx.ext.mathjax',  # LaTeX style math
    'sphinx.ext.viewcode',  # view code links
    'sphinx.ext.napoleon',  # for NumPy style docstrings
    'sphinx.ext.intersphinx',  # external links
    'sphinx.ext.autosummary',  # autosummary directive
    'sphinxext.custom_roles',  # local extension
    'sphinx_copybutton',
    'sphinx_automodapi.automodapi',  # see: https://github.com/lukelbd/sphinx-automodapi/tree/proplot-mods # noqa
    'nbsphinx',
]

extlinks = {
    'issue': ('https://github.com/lukelbd/proplot/issues/%s', 'GH#'),
    'commit': ('https://github.com/lukelbd/proplot/commit/%s', '@'),
    'pr': ('https://github.com/lukelbd/proplot/pull/%s', 'GH#'),
}

# Cupybutton configuration
# See: https://sphinx-copybutton.readthedocs.io/en/latest/
copybutton_prompt_text = r'>>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: '
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
copybutton_remove_prompts = True

# Fail on error. Note nbsphinx compiles *all* notebooks in docs unless excluded
nbsphinx_allow_errors = False

# Give *lots* of time for cell execution!
nbsphinx_timeout = 300

# Make nbsphinx detect jupytext
nbsphinx_custom_formats = {
    '.py': ['jupytext.reads', {'fmt': 'py:percent'}],
}

# Set InlineBackend params in case nbsphinx skips ones in config.py
# Not necessary because config.py configures the backend.
# nbsphinx_execute_arguments = [
#     "--InlineBackend.figure_formats={'svg'}",
#     "--InlineBackend.rc={'figure.dpi': 100}",
# ]

# For now do not run doctest tests, these are just to show syntax
# and expected output may be graphical
doctest_test_doctest_blocks = ''

# Generate stub pages whenever ::autosummary directive encountered
# This way don't have to call sphinx-autogen manually
autosummary_generate = True

# Use automodapi tool, created by astropy people. See:
# https://sphinx-automodapi.readthedocs.io/en/latest/automodapi.html#overview
# Normally have to *enumerate* function names manually. This will document
# them automatically. Just be careful, if you use from x import *, to exclude
# them in the automodapi:: directive
automodapi_toctreedirnm = 'api'  # create much better URL for the page
automodsumm_inherited_members = False

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
    'matplotlib': ('https://matplotlib.org/stable', None),
    'sphinx': ('http://www.sphinx-doc.org/en/stable', None),
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'xarray': ('http://xarray.pydata.org/en/stable', None),
    'cartopy': ('https://scitools.org.uk/cartopy/docs/latest', None),
    'basemap': ('https://matplotlib.org/basemap', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable', None),
    'pint': ('https://pint.readthedocs.io/en/stable/', None),
}

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False  # proplot imports everything in top-level namespace

# Napoleon options
# See: http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
# * use_param is set to False so that we can put multiple "parameters"
#   on one line -- for example 'xlocator, ylocator : locator-spec, optional'
# * docs claim napoleon_preprocess_types and napoleon_type_aliases only work
#   when napoleon_use_param is True but xarray sets to False and it still works
# * use_keyword is set to False because we do not want separate 'Keyword Arguments'
#   section and have same issue for multiple keywords.
# * use_ivar and use_rtype are set to False for (presumably) style consistency
#   with the above options set to False.
napoleon_use_ivar = False
napoleon_use_param = False
napoleon_use_keyword = False
napoleon_use_rtype = False
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False  # move init doc to 'class' doc
napoleon_preprocess_types = True
napoleon_type_aliases = {
    # General terms
    'sequence': ':term:`sequence`',
    'iterable': ':term:`iterable`',
    'callable': ':py:func:`callable`',
    'mapping': ':term:`mapping`',
    'hashable': ':term:`hashable <name>`',
    'dict-like': ':term:`dict-like <mapping>`',
    'path-like': ':term:`path-like <path-like object>`',
    'file-like': ':term:`file-like <file-like object>`',
    'bool': ':class:`bool <bool>`',
    'string': ':class:`string <str>`',
    'scalar': ':term:`scalar`',
    'array-like': ':term:`array-like`',
    # ProPlot terms
    'locator-like': ':term:`locator-like <locator>`',
    'formatter-like': ':term:`formatter-like <formatter>`',
    'scale-like': ':term:`scale-like <formatter>`',
    'colormap-like': ':term:`colormap-spec <colormap>`',
    'cycle-like': ':term:`cycle-spec <cycle>`',
    'color-like': ':py:func:`color-like <matplotlib.colors.is_color_like>`',
}

# Fix duplicate class member documentation from autosummary + numpydoc
# See: https://github.com/phn/pytpm/issues/3#issuecomment-12133978
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string.
source_suffix = ['.rst', '.html']

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
exclude_patterns = [
    '.DS_Store',
    '_build', '_templates', '_themes', 'trash', 'tmp',
    'conf.py', 'sphinxext', '*.ipynb', '**.ipynb_checkpoints',
]

# The name of the Pygments (syntax highlighting) style to use.
# The light-dark theme toggler overloads this, but set default anyway
pygments_style = 'none'

# Create local pygments copies
# Previously used: https://github.com/richleland/pygments-css
# But do not want to depend on some random repository
# from pygments.styles import get_all_styles  # noqa: E402
from pygments.formatters import HtmlFormatter  # noqa: E402
path = os.path.join('_static', 'pygments')
if not os.path.isdir(path):
    os.mkdir(path)
# for style in get_all_styles():
for style in ('pastie', 'monokai'):  # WARNING: update when _static/custom.js changes
    path = os.path.join('_static', 'pygments', style + '.css')
    if os.path.isfile(path):
        continue
    with open(path, 'w') as f:
        f.write(HtmlFormatter(style=style).get_style_defs('.highlight'))

# Create RST table and sample proplotrc file
from proplot.config import rc
rc._save_rst(os.path.join('_static', 'rctable.rst'))
rc._save_yaml(os.path.join('_static', 'proplotrc'))

# Role
# default family is py, but can also set default role so don't need
# :func:`name`, :module:`name`, etc.
default_role = 'py:obj'

# -- Options for HTML output -------------------------------------------------

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
# html_sidebars = {}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large. Static folder is for CSS and image files. Use ImageMagick to
# convert png to ico on command line with 'convert image.png image.ico'
html_favicon = 'logo_blank.ico'

# Logo
html_logo = 'logo_square.png'

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# Use modified RTD theme with overrides in custom.css and custom.js
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    'prev_next_buttons_location': 'bottom',  # top and bottom
}

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'proplotdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',

    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'proplot.tex', 'ProPlot Documentation',
     'Luke L. B. Davis', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (
        master_doc,
        'proplot',
        'ProPlot Documentation',
        [author],
        1
    )
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        'proplot',
        'ProPlot Documentation',
        author,
        'proplot',
        'A powerful matplotlib wrapper for making beautiful, '
        'publication-quality graphics.',
        'Miscellaneous'
    )
]


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True
