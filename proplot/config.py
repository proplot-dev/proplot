#!/usr/bin/env python3
"""
Tools for setting up proplot and configuring global settings.
See the :ref:`configuration guide <ug_config>` for details.
"""
# NOTE: The matplotlib analogue to this file is actually __init__.py
# but it makes more sense to have all the setup actions in a separate file
# so the namespace of the top-level module is unpolluted.
# NOTE: Why also load colormaps and cycles in this file and not colors.py?
# Because I think it makes sense to have all the code that "runs" (i.e. not
# just definitions) in the same place, and I was having issues with circular
# dependencies and where import order of __init__.py was affecting behavior.
import logging
import os
import re
import sys
from collections import namedtuple
from collections.abc import MutableMapping
from numbers import Real

import cycler
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
import matplotlib.mathtext  # noqa: F401
import matplotlib.style.core as mstyle
import numpy as np
from matplotlib import RcParams

from .internals import ic  # noqa: F401
from .internals import (
    _not_none,
    _pop_kwargs,
    _pop_props,
    _translate_grid,
    _version_mpl,
    docstring,
    rcsetup,
    warnings,
)

try:
    from IPython import get_ipython
except ImportError:

    def get_ipython():
        return


# Suppress warnings emitted by mathtext.py (_mathtext.py in recent versions)
# when when substituting dummy unavailable glyph due to fallback disabled.
logging.getLogger("matplotlib.mathtext").setLevel(logging.ERROR)

__all__ = [
    "Configurator",
    "rc",
    "rc_proplot",
    "rc_matplotlib",
    "use_style",
    "config_inline_backend",
    "register_cmaps",
    "register_cycles",
    "register_colors",
    "register_fonts",
    "RcConfigurator",  # deprecated
    "inline_backend_fmt",  # deprecated
]

# Constants
COLORS_KEEP = ("red", "green", "blue", "cyan", "yellow", "magenta", "white", "black")

# Configurator docstrings
_rc_docstring = """
local : bool, default: True
    Whether to load settings from the `~Configurator.local_files` file.
user : bool, default: True
    Whether to load settings from the `~Configurator.user_file` file.
default : bool, default: True
    Whether to reload built-in default proplot settings.
"""
docstring._snippet_manager["rc.params"] = _rc_docstring

# Registration docstrings
_shared_docstring = """
*args : path-spec or `~proplot.colors.{type}Colormap`, optional
    The {objects} to register. These can be file paths containing
    RGB data or `~proplot.colors.{type}Colormap` instances. By default,
    if positional arguments are passed, then `user` is set to ``False``.

    Valid file extensions are listed in the below table. Note that {objects}
    are registered according to their filenames -- for example, ``name.xyz``
    will be registered as ``'name'``.
"""  # noqa: E501
_cmap_exts_docstring = """
    ===================  ==========================================
    Extension            Description
    ===================  ==========================================
    ``.json``            JSON database of the channel segment data.
    ``.hex``             Comma-delimited list of HEX strings.
    ``.rgb``, ``.txt``   3-4 column table of channel values.
    ===================  ==========================================
"""
_cycle_exts_docstring = """
    ==================  ==========================================
    Extension           Description
    ==================  ==========================================
    ``.hex``            Comma-delimited list of HEX strings.
    ``.rgb``, ``.txt``  3-4 column table of channel values.
    ==================  ==========================================
"""
_color_docstring = """
*args : path-like or dict, optional
    The colors to register. These can be file paths containing RGB data or
    dictionary mappings of names to RGB values. By default, if positional
    arguments are passed, then `user` is set to ``False``. Files must have
    the extension ``.txt`` and should contain one line per color in the
    format ``name : hex``. Whitespace is ignored.
"""
_font_docstring = """
*args : path-like, optional
    The font files to add. By default, if positional arguments are passed, then
    `user` is set to ``False``. Files must have the extensions ``.ttf`` or ``.otf``.
    See `this link \
<https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting other font files to ``.ttf`` and ``.otf``.
"""
_register_docstring = """
user : bool, optional
    Whether to reload {objects} from `~Configurator.user_folder`. Default is
    ``False`` if positional arguments were passed and ``True`` otherwise.
local : bool, optional
    Whether to reload {objects} from `~Configurator.local_folders`. Default is
    ``False`` if positional arguments were passed and ``True`` otherwise.
default : bool, default: False
    Whether to reload the default {objects} packaged with proplot.
    Default is always ``False``.
"""
docstring._snippet_manager["rc.cmap_params"] = _register_docstring.format(
    objects="colormaps"
)  # noqa: E501
docstring._snippet_manager["rc.cycle_params"] = _register_docstring.format(
    objects="color cycles"
)  # noqa: E501
docstring._snippet_manager["rc.color_params"] = _register_docstring.format(
    objects="colors"
)  # noqa: E501
docstring._snippet_manager["rc.font_params"] = _register_docstring.format(
    objects="fonts"
)  # noqa: E501
docstring._snippet_manager["rc.cmap_args"] = _shared_docstring.format(
    objects="colormaps", type="Continuous"
)  # noqa: E501
docstring._snippet_manager["rc.cycle_args"] = _shared_docstring.format(
    objects="color cycles", type="Discrete"
)  # noqa: E501
docstring._snippet_manager["rc.color_args"] = _color_docstring
docstring._snippet_manager["rc.font_args"] = _font_docstring
docstring._snippet_manager["rc.cmap_exts"] = _cmap_exts_docstring
docstring._snippet_manager["rc.cycle_exts"] = _cycle_exts_docstring


def _init_user_file():
    """
    Initialize .proplotrc file.
    """
    file = Configurator.user_file()
    if not os.path.exists(file):
        Configurator._save_yaml(file, comment=True)


def _init_user_folders():
    """
    Initialize .proplot folder.
    """
    for subfolder in ("", "cmaps", "cycles", "colors", "fonts"):
        folder = Configurator.user_folder(subfolder)
        if not os.path.isdir(folder):
            os.mkdir(folder)


def _get_data_folders(folder, user=True, local=True, default=True, reverse=False):
    """
    Return data folder paths in reverse order of precedence.
    """
    # When loading colormaps, cycles, and colors, files in the latter
    # directories overwrite files in the former directories. When loading
    # fonts, the resulting paths need to be *reversed*.
    paths = []
    if default:
        paths.append(os.path.join(os.path.dirname(__file__), folder))
    if user:
        paths.append(Configurator.user_folder(folder))
    if local:
        paths.extend(Configurator.local_folders(folder))
    if reverse:
        paths = paths[::-1]
    return paths


def _iter_data_objects(folder, *args, **kwargs):
    """
    Iterate over input objects and files in the data folders that should be
    registered. Also yield an index indicating whether these are user files.
    """
    i = 0  # default files
    for i, path in enumerate(_get_data_folders(folder, **kwargs)):
        for dirname, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename[0] == ".":  # UNIX-style hidden files
                    continue
                path = os.path.join(dirname, filename)
                yield i, path
    i += 1  # user files
    for path in args:
        path = os.path.expanduser(path)
        if os.path.isfile(path):
            yield i, path
        else:
            raise FileNotFoundError(f"Invalid file path {path!r}.")


def _filter_style_dict(rcdict, warn=True):
    """
    Filter out blacklisted style parameters.
    """
    # NOTE: This implements bugfix: https://github.com/matplotlib/matplotlib/pull/17252
    # This fix is *critical* for proplot because we always run style.use()
    # when the configurator is made. Without fix backend is reset every time
    # you import proplot in jupyter notebooks. So apply retroactively.
    rcdict_filtered = {}
    for key in rcdict:
        if key in mstyle.STYLE_BLACKLIST:
            if warn:
                warnings._warn_proplot(
                    f"Dictionary includes a parameter, {key!r}, that is not related "
                    "to style. Ignoring."
                )
        else:
            rcdict_filtered[key] = rcdict[key]
    return rcdict_filtered


def _get_default_style_dict():
    """
    Get the default rc parameters dictionary with deprecated parameters filtered.
    """
    # NOTE: Use RcParams update to filter and translate deprecated settings
    # before actually applying them to rcParams down pipeline. This way we can
    # suppress warnings for deprecated default params but still issue warnings
    # when user-supplied stylesheets have deprecated params.
    # WARNING: Some deprecated rc params remain in dictionary as None so we
    # filter them out. Beware if hidden attribute changes.
    # WARNING: The examples.directory deprecation was handled specially inside
    # RcParams in early versions. Manually pop it out here.
    rcdict = _filter_style_dict(mpl.rcParamsDefault, warn=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", mpl.MatplotlibDeprecationWarning)
        rcdict = dict(RcParams(rcdict))
    for attr in ("_deprecated_set", "_deprecated_remain_as_none"):
        deprecated = getattr(mpl, attr, ())
        for key in deprecated:  # _deprecated_set is in matplotlib < 3.4
            rcdict.pop(key, None)
    rcdict.pop("examples.directory", None)  # special case for matplotlib < 3.2
    return rcdict


def _get_style_dict(style, filter=True):
    """
    Return a dictionary of settings belonging to the requested style(s). If `filter`
    is ``True``, invalid style parameters like `backend` are filtered out.
    """
    # NOTE: This is adapted from matplotlib source for the following changes:
    # 1. Add an 'original' pseudo style. Like rcParamsOrig except we also reload
    #    from the user matplotlibrc file.
    # 2. When the style is changed we reset to the default state ignoring matplotlibrc.
    #    By contrast matplotlib applies styles on top of current state (including
    #    matplotlibrc changes and runtime rcParams changes) but the word 'style'
    #    implies a rigid static format. This makes more sense.
    # 3. Add a separate function that returns lists of style dictionaries so that
    #    we can modify the active style in a context block. Proplot context is more
    #    conservative than matplotlib's rc_context because it gets called a lot
    #    (e.g. every time you make an axes and every format() call). Instead of
    #    copying the entire rcParams dict we just track the keys that were changed.
    style_aliases = {
        "538": "fivethirtyeight",
        "mpl20": "default",
        "mpl15": "classic",
        "original": mpl.matplotlib_fname(),
        "seaborn": "seaborn-v0_8",
    }

    # Always apply the default style *first* so styles are rigid
    kw_matplotlib = _get_default_style_dict()
    if style == "default" or style is mpl.rcParamsDefault:
        return kw_matplotlib

    # Apply limited deviations from the matplotlib style that we want to propagate to
    # other styles. Want users selecting stylesheets to have few surprises, so
    # currently just enforce the new aesthetically pleasing fonts.
    kw_matplotlib["font.family"] = "sans-serif"
    for fmly in ("serif", "sans-serif", "monospace", "cursive", "fantasy"):
        kw_matplotlib["font." + fmly] = rcsetup._rc_matplotlib_default["font." + fmly]

    # Apply user input style(s) one by one
    if isinstance(style, str) or isinstance(style, dict):
        styles = [style]
    else:
        styles = style
    for style in styles:
        if isinstance(style, dict):
            kw = style
        elif isinstance(style, str):
            style = style_aliases.get(style, style)
            if style in mstyle.library:
                kw = mstyle.library[style]
            else:
                try:
                    kw = mpl.rc_params_from_file(style, use_default_template=False)
                except IOError:
                    raise IOError(
                        f"Style {style!r} not found in the style library and input "
                        "is not a valid URL or file path. Available styles are: "
                        + ", ".join(map(repr, mstyle.available))
                        + "."
                    )
        else:
            raise ValueError(f"Invalid style {style!r}. Must be string or dictionary.")
        if filter:
            kw = _filter_style_dict(kw, warn=True)
        kw_matplotlib.update(kw)

    return kw_matplotlib


def _infer_proplot_dict(kw_params):
    """
    Infer values for proplot's "added" parameters from stylesheet parameters.
    """
    kw_proplot = {}
    mpl_to_proplot = {
        "xtick.labelsize": (
            "tick.labelsize",
            "grid.labelsize",
        ),
        "ytick.labelsize": (
            "tick.labelsize",
            "grid.labelsize",
        ),
        "axes.titlesize": (
            "abc.size",
            "suptitle.size",
            "title.size",
            "leftlabel.size",
            "rightlabel.size",
            "toplabel.size",
            "bottomlabel.size",
        ),
        "text.color": (
            "abc.color",
            "suptitle.color",
            "title.color",
            "tick.labelcolor",
            "grid.labelcolor",
            "leftlabel.color",
            "rightlabel.color",
            "toplabel.color",
            "bottomlabel.color",
        ),
    }
    for key, params in mpl_to_proplot.items():
        if key in kw_params:
            value = kw_params[key]
            for param in params:
                kw_proplot[param] = value
    return kw_proplot


def config_inline_backend(fmt=None):
    """
    Set up the ipython `inline backend display format \
<https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    and ensure that inline figures always look the same as saved figures.
    This runs the following ipython magic commands:

    .. code-block:: ipython

        %config InlineBackend.figure_formats = rc['inlineformat']
        %config InlineBackend.rc = {}  # never override rc settings
        %config InlineBackend.close_figures = True  # cells start with no active figures
        %config InlineBackend.print_figure_kwargs = {'bbox_inches': None}

    When the inline backend is inactive or unavailable, this has no effect.
    This function is called when you modify the :rcraw:`inlineformat` property.

    Parameters
    ----------
    fmt : str or sequence, default: :rc:`inlineformat`
        The inline backend file format or a list thereof. Valid formats
        include ``'jpg'``, ``'png'``, ``'svg'``, ``'pdf'``, and ``'retina'``.

    See also
    --------
    Configurator
    """
    # Note if inline backend is unavailable this will fail silently
    ipython = get_ipython()
    if ipython is None:
        return
    fmt = _not_none(fmt, rc_proplot["inlineformat"])
    if isinstance(fmt, str):
        fmt = [fmt]
    elif np.iterable(fmt):
        fmt = list(fmt)
    else:
        raise ValueError(f"Invalid inline backend format {fmt!r}. Must be string.")
    ipython.magic("config InlineBackend.figure_formats = " + repr(fmt))
    ipython.magic("config InlineBackend.rc = {}")
    ipython.magic("config InlineBackend.close_figures = True")
    ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches': None}")


def use_style(style):
    """
    Apply the `matplotlib style(s) \
<https://matplotlib.org/stable/tutorials/introductory/customizing.html>`__
    with `matplotlib.style.use`. This function is
    called when you modify the :rcraw:`style` property.

    Parameters
    ----------
    style : str or sequence or dict-like
        The matplotlib style name(s) or stylesheet filename(s), or dictionary(s)
        of settings. Use ``'default'`` to apply matplotlib default settings and
        ``'original'`` to include settings from your user ``matplotlibrc`` file.

    See also
    --------
    Configurator
    matplotlib.style.use
    """
    # NOTE: This function is not really necessary but makes proplot's
    # stylesheet-supporting features obvious. Plus changing the style does
    # so much *more* than changing rc params or quick settings, so it is
    # nice to have dedicated function instead of just another rc_param name.
    kw_matplotlib = _get_style_dict(style)
    rc_matplotlib.update(kw_matplotlib)
    rc_proplot.update(_infer_proplot_dict(kw_matplotlib))


@docstring._snippet_manager
def register_cmaps(*args, user=None, local=None, default=False):
    """
    Register named colormaps. This is called on import.

    Parameters
    ----------
    %(rc.cmap_args)s

    %(rc.cmap_exts)s

    %(rc.cmap_params)s

    See also
    --------
    register_cycles
    register_colors
    register_fonts
    proplot.demos.show_cmaps
    """
    # Register input colormaps
    from . import colors as pcolors

    user = _not_none(user, not bool(args))  # skip user folder if input args passed
    local = _not_none(local, not bool(args))
    paths = []
    for arg in args:
        if isinstance(arg, mcolors.Colormap):
            pcolors._cmap_database.register(arg, name=arg.name)
        else:
            paths.append(arg)

    # Register data files
    for i, path in _iter_data_objects(
        "cmaps", *paths, user=user, local=local, default=default
    ):
        cmap = pcolors.ContinuousColormap.from_file(path, warn_on_failure=True)
        if not cmap:
            continue
        if i == 0 and cmap.name.lower() in pcolors.CMAPS_CYCLIC:
            cmap.set_cyclic(True)
        pcolors._cmap_database.register(cmap, name=cmap.name)


@docstring._snippet_manager
def register_cycles(*args, user=None, local=None, default=False):
    """
    Register named color cycles. This is called on import.

    Parameters
    ----------
    %(rc.cycle_args)s

    %(rc.cycle_exts)s

    %(rc.cycle_params)s

    See also
    --------
    register_cmaps
    register_colors
    register_fonts
    proplot.demos.show_cycles
    """
    # Register input color cycles
    from . import colors as pcolors

    user = _not_none(user, not bool(args))  # skip user folder if input args passed
    local = _not_none(local, not bool(args))
    paths = []
    for arg in args:
        if isinstance(arg, mcolors.Colormap):
            pcolors._cmap_database.register(arg, name=arg.name)
        else:
            paths.append(arg)

    # Register data files
    for _, path in _iter_data_objects(
        "cycles", *paths, user=user, local=local, default=default
    ):
        cmap = pcolors.DiscreteColormap.from_file(path, warn_on_failure=True)
        if not cmap:
            continue
        pcolors._cmap_database.register(cmap, name=cmap.name)


@docstring._snippet_manager
def register_colors(
    *args, user=None, local=None, default=False, space=None, margin=None, **kwargs
):
    """
    Register named colors. This is called on import.

    Parameters
    ----------
    %(rc.color_args)s
    %(rc.color_params)s
    space : {'hcl', 'hsl', 'hpl'}, optional
        The colorspace used to pick "perceptually distinct" colors from
        the `XKCD color survey <https://xkcd.com/color/rgb/>`__.
        If passed then `default` is set to ``True``.
    margin : float, default: 0.1
        The margin used to pick "perceptually distinct" colors from the
        `XKCD color survey <https://xkcd.com/color/rgb/>`__. The normalized hue,
        saturation, and luminance of each color must differ from the channel
        values of the prededing colors by `margin` in order to be registered.
        Must fall between ``0`` and ``1`` (``0`` will register all colors).
        If passed then `default` is set to ``True``.
    **kwargs
        Additional color name specifications passed as keyword arguments rather
        than positional argument dictionaries.

    See also
    --------
    register_cmaps
    register_cycles
    register_fonts
    proplot.demos.show_colors
    """
    from . import colors as pcolors

    default = default or space is not None or margin is not None
    margin = _not_none(margin, 0.1)
    space = _not_none(space, "hcl")

    # Remove previously registered colors
    # NOTE: Try not to touch matplotlib colors for compatibility
    srcs = {"xkcd": pcolors.COLORS_XKCD, "opencolor": pcolors.COLORS_OPEN}
    if default:  # possibly slow but not these dicts are empty on startup
        for src in srcs.values():
            for key in src:
                if key not in COLORS_KEEP:
                    pcolors._color_database.pop(key, None)  # this also clears cache
            src.clear()

    # Register input colors
    user = _not_none(user, not bool(args) and not bool(kwargs))  # skip if args passed
    local = _not_none(local, not bool(args) and not bool(kwargs))
    paths = []
    for arg in args:
        if isinstance(arg, dict):
            kwargs.update(arg)
        else:
            paths.append(arg)
    for key, color in kwargs.items():
        if mcolors.is_color_like(color):
            pcolors._color_database[key] = mcolors.to_rgba(color)
        else:
            raise ValueError(f"Invalid color {key}={color!r}.")

    # Load colors from file and get their HCL values
    # NOTE: Colors that come *later* overwrite colors that come earlier.
    for i, path in _iter_data_objects(
        "colors", *paths, user=user, local=local, default=default
    ):
        loaded = pcolors._load_colors(path, warn_on_failure=True)
        if i == 0:
            cat, _ = os.path.splitext(os.path.basename(path))
            if cat not in srcs:
                raise RuntimeError(f"Unknown proplot color database {path!r}.")
            src = srcs[cat]
            if cat == "xkcd":
                for key in COLORS_KEEP:
                    loaded[key] = pcolors._color_database[key]  # keep the same
                loaded = pcolors._standardize_colors(loaded, space, margin)
            src.clear()
            src.update(loaded)  # needed for demos.show_colors()
        pcolors._color_database.update(loaded)


@docstring._snippet_manager
def register_fonts(*args, user=True, local=True, default=False):
    """
    Register font families. This is called on import.

    Parameters
    ----------
    %(rc.font_args)s
    %(rc.font_params)s

    See also
    --------
    register_cmaps
    register_cycles
    register_colors
    proplot.demos.show_fonts
    """
    # Find proplot fonts
    # WARNING: Must search data files in reverse because font manager will
    # not overwrite existing fonts with user-input fonts.
    # WARNING: If you include a font file with an unrecognized style,
    # matplotlib may use that font instead of the 'normal' one! Valid styles:
    # 'ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman',
    # 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'
    # https://matplotlib.org/api/font_manager_api.html
    # For macOS the only fonts with 'Thin' in one of the .ttf file names
    # are Helvetica Neue and .SF NS Display Condensed. Never try to use these!
    paths_proplot = _get_data_folders(
        "fonts", user=user, local=local, default=default, reverse=True
    )
    fnames_proplot = set(mfonts.findSystemFonts(paths_proplot))
    for path in args:
        path = os.path.expanduser(path)
        if os.path.isfile(path):
            fnames_proplot.add(path)
        else:
            raise FileNotFoundError(f"Invalid font file path {path!r}.")

    # Detect user-input ttc fonts and issue warning
    fnames_proplot_ttc = {
        file for file in fnames_proplot if os.path.splitext(file)[1] == ".ttc"
    }
    if fnames_proplot_ttc:
        warnings._warn_proplot(
            "Ignoring the following .ttc fonts because they cannot be "
            "saved into PDF or EPS files (see matplotlib issue #3135): "
            + ", ".join(map(repr, sorted(fnames_proplot_ttc)))
            + ". Please consider expanding them into separate .ttf files."
        )

    # Rebuild font cache only if necessary! Can be >50% of total import time!
    fnames_all = {font.fname for font in mfonts.fontManager.ttflist}
    fnames_proplot -= fnames_proplot_ttc
    if not fnames_all >= fnames_proplot:
        warnings._warn_proplot(
            "Rebuilding font cache. This usually happens "
            "after installing or updating proplot."
        )
        if hasattr(mfonts.fontManager, "addfont"):
            # Newer API lets us add font files manually and deprecates TTFPATH. However
            # to cache fonts added this way, we must call json_dump explicitly.
            # NOTE: Previously, cache filename was specified as _fmcache variable, but
            # recently became inaccessible. Must reproduce mpl code instead.
            # NOTE: Older mpl versions used fontList.json as the cache, but these
            # versions also did not have 'addfont', so makes no difference.
            for fname in fnames_proplot:
                mfonts.fontManager.addfont(fname)
            cache = os.path.join(
                mpl.get_cachedir(), f"fontlist-v{mfonts.FontManager.__version__}.json"
            )
            mfonts.json_dump(mfonts.fontManager, cache)
        else:
            # Older API requires us to modify TTFPATH
            # NOTE: Previously we tried to modify TTFPATH before importing
            # font manager with hope that it would load proplot fonts on
            # initialization. But 99% of the time font manager just imports
            # the FontManager from cache, so we would have to rebuild anyway.
            paths = ":".join(paths_proplot)
            if "TTFPATH" not in os.environ:
                os.environ["TTFPATH"] = paths
            elif paths not in os.environ["TTFPATH"]:
                os.environ["TTFPATH"] += ":" + paths
            mfonts._rebuild()

    # Remove ttc files and 'Thin' fonts *after* rebuild
    # NOTE: 'Thin' filter is ugly kludge but without this matplotlib picks up on
    # Roboto thin ttf files installed on the RTD server when compiling docs.
    mfonts.fontManager.ttflist = [
        font
        for font in mfonts.fontManager.ttflist
        if os.path.splitext(font.fname)[1] != ".ttc"
        and (_version_mpl >= "3.3" or "Thin" not in os.path.basename(font.fname))
    ]


class Configurator(MutableMapping, dict):
    """
    A dictionary-like class for managing `matplotlib settings
    <https://matplotlib.org/stable/tutorials/introductory/customizing.html>`__
    stored in `rc_matplotlib` and :ref:`proplot settings <ug_rcproplot>`
    stored in `rc_proplot`. This class is instantiated as the `rc` object
    on import. See the :ref:`user guide <ug_config>` for details.
    """

    def __repr__(self):
        cls = type("rc", (dict,), {})  # temporary class with short name
        src = cls({key: val for key, val in rc_proplot.items() if "." not in key})
        return type(rc_matplotlib).__repr__(src).strip()[:-1] + ",\n    ...\n    })"

    def __str__(self):
        cls = type("rc", (dict,), {})  # temporary class with short name
        src = cls({key: val for key, val in rc_proplot.items() if "." not in key})
        return type(rc_matplotlib).__str__(src) + "\n..."

    def __iter__(self):
        yield from rc_proplot  # sorted proplot settings, ignoring deprecations
        yield from rc_matplotlib  # sorted matplotlib settings, ignoring deprecations

    def __len__(self):
        return len(tuple(iter(self)))

    def __delitem__(self, key):  # noqa: U100
        raise RuntimeError("rc settings cannot be deleted.")

    def __delattr__(self, attr):  # noqa: U100
        raise RuntimeError("rc settings cannot be deleted.")

    @docstring._snippet_manager
    def __init__(self, local=True, user=True, default=True, **kwargs):
        """
        Parameters
        ----------
        %(rc.params)s
        """
        self._context = []
        self._init(local=local, user=user, default=default, **kwargs)

    def __getitem__(self, key):
        """
        Return an `rc_matplotlib` or `rc_proplot` setting using dictionary notation
        (e.g., ``value = pplt.rc[name]``).
        """
        key, _ = self._validate_key(key)  # might issue proplot removed/renamed error
        try:
            return rc_proplot[key]
        except KeyError:
            pass
        return rc_matplotlib[key]  # might issue matplotlib removed/renamed error

    def __setitem__(self, key, value):
        """
        Modify an `rc_matplotlib` or `rc_proplot` setting using dictionary notation
        (e.g., ``pplt.rc[name] = value``).
        """
        kw_proplot, kw_matplotlib = self._get_item_dicts(key, value)
        rc_proplot.update(kw_proplot)
        rc_matplotlib.update(kw_matplotlib)

    def __getattr__(self, attr):
        """
        Return an `rc_matplotlib` or `rc_proplot` setting using "dot" notation
        (e.g., ``value = pplt.rc.name``).
        """
        if attr[:1] == "_":
            return super().__getattribute__(attr)  # raise built-in error
        else:
            return self.__getitem__(attr)

    def __setattr__(self, attr, value):
        """
        Modify an `rc_matplotlib` or `rc_proplot` setting using "dot" notation
        (e.g., ``pplt.rc.name = value``).
        """
        if attr[:1] == "_":
            super().__setattr__(attr, value)
        else:
            self.__setitem__(attr, value)

    def __enter__(self):
        """
        Apply settings from the most recent context block.
        """
        if not self._context:
            raise RuntimeError(
                "rc object must be initialized for context block using rc.context()."
            )
        context = self._context[-1]
        kwargs = context.kwargs
        rc_new = context.rc_new  # used for context-based _get_item_context
        rc_old = context.rc_old  # used to re-apply settings without copying whole dict
        for key, value in kwargs.items():
            kw_proplot, kw_matplotlib = self._get_item_dicts(key, value)
            for rc_dict, kw_new in zip(
                (rc_proplot, rc_matplotlib),
                (kw_proplot, kw_matplotlib),
            ):
                for key, value in kw_new.items():
                    rc_old[key] = rc_dict[key]
                    rc_new[key] = rc_dict[key] = value

    def __exit__(self, *args):  # noqa: U100
        """
        Restore settings from the most recent context block.
        """
        if not self._context:
            raise RuntimeError(
                "rc object must be initialized for context block using rc.context()."
            )
        context = self._context[-1]
        for key, value in context.rc_old.items():
            kw_proplot, kw_matplotlib = self._get_item_dicts(key, value)
            rc_proplot.update(kw_proplot)
            rc_matplotlib.update(kw_matplotlib)
        del self._context[-1]

    def _init(self, *, local, user, default, skip_cycle=False):
        """
        Initialize the configurator.
        """
        # Always remove context objects
        self._context.clear()

        # Update from default settings
        # NOTE: see _remove_blacklisted_style_params bugfix
        if default:
            rc_matplotlib.update(_get_style_dict("original", filter=False))
            rc_matplotlib.update(rcsetup._rc_matplotlib_default)
            rc_proplot.update(rcsetup._rc_proplot_default)
            for key, value in rc_proplot.items():
                kw_proplot, kw_matplotlib = self._get_item_dicts(
                    key, value, skip_cycle=skip_cycle
                )
                rc_matplotlib.update(kw_matplotlib)
                rc_proplot.update(kw_proplot)

        # Update from user home
        user_path = None
        if user:
            user_path = self.user_file()
            if os.path.isfile(user_path):
                self.load(user_path)

        # Update from local paths
        if local:
            local_paths = self.local_files()
            for path in local_paths:
                if path == user_path:  # local files always have precedence
                    continue
                self.load(path)

    @staticmethod
    def _validate_key(key, value=None):
        """
        Validate setting names and handle `rc_proplot` deprecations.
        """
        # NOTE: Not necessary to check matplotlib key here because... not sure why.
        # Think deprecated matplotlib keys are not involved in any synced settings.
        # Also note _check_key includes special handling for some renamed keys.
        if not isinstance(key, str):
            raise KeyError(f"Invalid key {key!r}. Must be string.")
        key = key.lower()
        if "." not in key:
            key = rcsetup._rc_nodots.get(key, key)
        key, value = rc_proplot._check_key(key, value)  # may issue deprecation warning
        return key, value

    @staticmethod
    def _validate_value(key, value):
        """
        Validate setting values and convert numpy ndarray to list if possible.
        """
        # NOTE: Ideally would implicitly validate on subsequent assignment to rc
        # dictionaries, but must explicitly do it here, so _get_item_dicts can
        # work with e.g. 'tick.lenratio', so _get_item_dicts does not have to include
        # deprecated name handling in its if statements, and so _load_file can
        # catch errors and emit warnings with line number indications as files
        # are being read rather than after the end of the file reading.
        if isinstance(value, np.ndarray):
            value = value.item() if value.size == 1 else value.tolist()
        validate_matplotlib = getattr(rc_matplotlib, "validate", None)
        validate_proplot = rc_proplot._validate
        if validate_matplotlib is not None and key in validate_matplotlib:
            value = validate_matplotlib[key](value)
        elif key in validate_proplot:
            value = validate_proplot[key](value)
        return value

    def _get_item_context(self, key, mode=None):
        """
        As with `~Configurator.__getitem__` but the search is limited based
        on the context mode and ``None`` is returned if the key is not found.
        """
        key, _ = self._validate_key(key)
        if mode is None:
            mode = self._context_mode
        cache = tuple(context.rc_new for context in self._context)
        if mode == 0:
            rcdicts = (*cache, rc_proplot, rc_matplotlib)
        elif mode == 1:
            rcdicts = (*cache, rc_proplot)  # added settings only!
        elif mode == 2:
            rcdicts = (*cache,)
        else:
            raise ValueError(f"Invalid caching mode {mode!r}.")
        for rcdict in rcdicts:
            if not rcdict:
                continue
            try:
                return rcdict[key]
            except KeyError:
                continue
        if mode == 0:  # otherwise return None
            raise KeyError(f"Invalid rc setting {key!r}.")

    def _get_item_dicts(self, key, value, skip_cycle=False):
        """
        Return dictionaries for updating the `rc_proplot` and `rc_matplotlib`
        properties associated with this key. Used when setting items, entering
        context blocks, or loading files.
        """
        # Get validated key, value, and child keys
        key, value = self._validate_key(key, value)
        value = self._validate_value(key, value)
        keys = (key,) + rcsetup._rc_children.get(key, ())  # settings to change
        contains = lambda *args: any(arg in keys for arg in args)  # noqa: E731

        # Fill dictionaries of matplotlib and proplot settings
        # NOTE: Raise key error right away so it can be caught by _load_file().
        # Also ignore deprecation warnings so we only get them *once* on assignment
        kw_proplot = {}  # custom properties
        kw_matplotlib = {}  # builtin properties
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", mpl.MatplotlibDeprecationWarning)
            warnings.simplefilter("ignore", warnings.ProplotWarning)
            for key in keys:
                if key in rc_matplotlib:
                    kw_matplotlib[key] = value
                elif key in rc_proplot:
                    kw_proplot[key] = value
                else:
                    raise KeyError(f"Invalid rc setting {key!r}.")

        # Special key: configure inline backend
        if contains("inlineformat"):
            config_inline_backend(value)

        # Special key: apply stylesheet
        elif contains("style"):
            if value is not None:
                ikw_matplotlib = _get_style_dict(value)
                kw_matplotlib.update(ikw_matplotlib)
                kw_proplot.update(_infer_proplot_dict(ikw_matplotlib))

        # Cycler
        # NOTE: Have to skip this step during initial proplot import
        elif contains("cycle") and not skip_cycle:
            from .colors import _get_cmap_subtype

            cmap = _get_cmap_subtype(value, "discrete")
            kw_matplotlib["axes.prop_cycle"] = cycler.cycler("color", cmap.colors)
            kw_matplotlib["patch.facecolor"] = "C0"

        # Turning bounding box on should turn border off and vice versa
        elif contains("abc.bbox", "title.bbox", "abc.border", "title.border"):
            if value:
                name, this = key.split(".")
                other = "border" if this == "bbox" else "bbox"
                kw_proplot[name + "." + other] = False

        # Fontsize
        # NOTE: Re-application of e.g. size='small' uses the updated 'font.size'
        elif contains("font.size"):
            kw_proplot.update(
                {
                    key: value
                    for key, value in rc_proplot.items()
                    if key in rcsetup.FONT_KEYS and value in mfonts.font_scalings
                }
            )
            kw_matplotlib.update(
                {
                    key: value
                    for key, value in rc_matplotlib.items()
                    if key in rcsetup.FONT_KEYS and value in mfonts.font_scalings
                }
            )

        # Tick length/major-minor tick length ratio
        elif contains("tick.len", "tick.lenratio"):
            if contains("tick.len"):
                ticklen = value
                ratio = rc_proplot["tick.lenratio"]
            else:
                ticklen = rc_proplot["tick.len"]
                ratio = value
            kw_matplotlib["xtick.minor.size"] = ticklen * ratio
            kw_matplotlib["ytick.minor.size"] = ticklen * ratio

        # Spine width/major-minor tick width ratio
        elif contains("tick.width", "tick.widthratio"):
            if contains("tick.width"):
                tickwidth = value
                ratio = rc_proplot["tick.widthratio"]
            else:
                tickwidth = rc_proplot["tick.width"]
                ratio = value
            kw_matplotlib["xtick.minor.width"] = tickwidth * ratio
            kw_matplotlib["ytick.minor.width"] = tickwidth * ratio

        # Gridline width
        elif contains("grid.width", "grid.widthratio"):
            if contains("grid.width"):
                gridwidth = value
                ratio = rc_proplot["grid.widthratio"]
            else:
                gridwidth = rc_proplot["grid.width"]
                ratio = value
            kw_proplot["gridminor.linewidth"] = gridwidth * ratio
            kw_proplot["gridminor.width"] = gridwidth * ratio

        # Gridline toggling
        elif contains("grid", "gridminor"):
            b, which = _translate_grid(
                value, "gridminor" if contains("gridminor") else "grid"
            )
            kw_matplotlib["axes.grid"] = b
            kw_matplotlib["axes.grid.which"] = which

        return kw_proplot, kw_matplotlib

    @staticmethod
    def _get_axisbelow_zorder(axisbelow):
        """
        Convert the `axisbelow` string to its corresponding `zorder`.
        """
        if axisbelow is True:
            zorder = 0.5
        elif axisbelow is False:
            zorder = 2.5
        elif axisbelow in ("line", "lines"):
            zorder = 1.5
        else:
            raise ValueError(f"Unexpected axisbelow value {axisbelow!r}.")
        return zorder

    def _get_background_props(self, patch_kw=None, native=True, **kwargs):
        """
        Return background properties, optionally filtering the output dictionary
        based on the context.
        """
        # Deprecated behavior
        context = native or self._context_mode == 2
        if patch_kw:
            warnings._warn_proplot(
                "'patch_kw' is no longer necessary as of proplot v0.8. "
                "Pass the parameters as keyword arguments instead."
            )
            kwargs.update(patch_kw)

        # Get user-input properties and changed rc settings
        # NOTE: Here we use 'color' as an alias for just 'edgecolor' rather than
        # both 'edgecolor' and 'facecolor' to match 'xcolor' and 'ycolor' arguments.
        props = _pop_props(kwargs, "patch")
        if "color" in props:
            props.setdefault("edgecolor", props.pop("color"))
        for key in ("alpha", "facecolor", "linewidth", "edgecolor"):
            value = self.find("axes." + key, context=context)
            if value is not None:
                props.setdefault(key, value)

        # Partition properties into face and edge
        kw_face = _pop_kwargs(props, "alpha", "facecolor")
        kw_edge = _pop_kwargs(props, "edgecolor", "linewidth", "linestyle")
        kw_edge["capstyle"] = "projecting"  # NOTE: needed to fix cartopy bounds
        if "color" in props:
            kw_edge.setdefault("edgecolor", props.pop("color"))
        if kwargs:
            raise TypeError(f"Unexpected keyword argument(s): {kwargs!r}")

        return kw_face, kw_edge

    def _get_gridline_bool(self, grid=None, axis=None, which="major", native=True):
        """
        Return major and minor gridline toggles from ``axes.grid``, ``axes.grid.which``,
        and ``axes.grid.axis``, optionally returning `None` based on the context.
        """
        # NOTE: If you pass 'grid' or 'gridminor' the native args are updated
        # NOTE: Very careful to return not None only if setting was changed.
        # Avoid unnecessarily triggering grid redraws (esp. bad for geo.py)
        context = native or self._context_mode == 2
        grid_on = self.find("axes.grid", context=context)
        which_on = self.find("axes.grid.which", context=context)
        if grid_on is not None or which_on is not None:  # if *one* was changed
            axis_on = self["axes.grid.axis"]  # always need this property
            grid_on = _not_none(grid_on, self["axes.grid"])
            which_on = _not_none(which_on, self["axes.grid.which"])
            axis = _not_none(axis, "x")
            axis_on = axis is None or axis_on in (axis, "both")
            which_on = which_on in (which, "both")
            grid = _not_none(grid, grid_on and axis_on and which_on)
        return grid

    def _get_gridline_props(self, which="major", native=True, rebuild=False):
        """
        Return gridline properties, optionally filtering the output dictionary
        based on the context.
        """
        # Line properties
        # NOTE: Gridline zorder is controlled automatically by matplotlib but
        # must be controlled manually for geographic projections
        key = "grid" if which == "major" else "gridminor"
        prefix = "grid_" if native else ""  # for native gridlines use this prefix
        context = not rebuild and (native or self._context_mode == 2)
        kwlines = self.fill(
            {
                f"{prefix}alpha": f"{key}.alpha",
                f"{prefix}color": f"{key}.color",
                f"{prefix}linewidth": f"{key}.linewidth",
                f"{prefix}linestyle": f"{key}.linestyle",
            },
            context=context,
        )
        axisbelow = self.find("axes.axisbelow", context=context)
        if axisbelow is not None:
            if native:  # this is a native plot so use set_axisbelow() down the line
                kwlines["axisbelow"] = axisbelow
            else:  # this is a geographic plot so apply with zorder
                kwlines["zorder"] = self._get_axisbelow_zorder(axisbelow)
        return kwlines

    def _get_label_props(self, native=True, **kwargs):
        """
        Return the axis label properties, optionally filtering the output dictionary
        based on the context.
        """
        # Get the label settings
        # NOTE: This permits passing arbitrary additional args to set_[xy]label()
        context = native or self._context_mode == 2
        kw = self.fill(
            {
                "color": "axes.labelcolor",
                "weight": "axes.labelweight",
                "size": "axes.labelsize",
                "family": "font.family",
                "labelpad": "axes.labelpad",  # read by set_xlabel/set_ylabel
            },
            context=context,
        )
        for key, value in kwargs.items():
            if value is not None:  # allow e.g. color=None
                kw[key] = value
        return kw

    def _get_loc_string(self, string, axis=None, native=True):
        """
        Return `tickloc` and `spineloc` location strings from the `rc` boolean toggles,
        optionally returning `None` based on the context.
        """
        context = native or self._context_mode == 2
        axis = _not_none(axis, "x")
        opt1, opt2 = ("top", "bottom") if axis == "x" else ("left", "right")
        b1 = self.find(f"{string}.{opt1}", context=context)
        b2 = self.find(f"{string}.{opt2}", context=context)
        if b1 is None and b2 is None:
            return None
        elif b1 and b2:
            return "both"
        elif b1:
            return opt1
        elif b2:
            return opt2
        else:
            return "neither"

    def _get_tickline_props(self, axis=None, which="major", native=True, rebuild=False):
        """
        Return the tick line properties, optionally filtering the output dictionary
        based on the context.
        """
        # Tick properties obtained with rc.category
        # NOTE: This loads 'size', 'width', 'pad', 'bottom', and 'top'
        axis = _not_none(axis, "x")
        context = not rebuild and (native or self._context_mode == 2)
        kwticks = self.category(f"{axis}tick.{which}", context=context)
        kwticks.pop("visible", None)
        for key in ("color", "direction"):
            value = self.find(f"{axis}tick.{key}", context=context)
            if value is not None:
                kwticks[key] = value
        return kwticks

    def _get_ticklabel_props(self, axis=None, native=True, rebuild=False):
        """
        Return the tick label properties, optionally filtering the output dictionary
        based on the context.
        """
        # NOTE: 'tick.label' properties are now synonyms of 'grid.label' properties
        sprefix = axis or ""
        cprefix = sprefix if _version_mpl >= "3.4" else ""  # new settings
        context = not rebuild and (native or self._context_mode == 2)
        kwtext = self.fill(
            {
                "color": f"{cprefix}tick.labelcolor",  # native setting sometimes avail
                "size": f"{sprefix}tick.labelsize",  # native setting always avail
                "weight": "tick.labelweight",  # native setting never avail
                "family": "font.family",  # apply manually
            },
            context=context,
        )
        if kwtext.get("color", None) == "inherit":
            # Inheritence is not automatic for geographic
            # gridline labels so we apply inheritence here.
            kwtext["color"] = self[f"{sprefix}tick.color"]
        return kwtext

    @staticmethod
    def local_files():
        """
        Return locations of files named ``proplotrc`` in this directory and in parent
        directories. "Hidden" files with a leading dot are also recognized. These are
        automatically loaded when proplot is imported.

        See also
        --------
        Configurator.user_file
        Configurator.local_folders
        """
        cdir = os.getcwd()
        paths = []
        while cdir:  # i.e. not root
            for name in ("proplotrc", ".proplotrc"):
                path = os.path.join(cdir, name)
                if os.path.isfile(path):
                    paths.append(path)
            ndir = os.path.dirname(cdir)
            if ndir == cdir:  # root
                break
            cdir = ndir
        return paths[::-1]  # sort from decreasing to increasing importantce

    @staticmethod
    def local_folders(subfolder=None):
        """
        Return locations of folders named ``proplot_cmaps``, ``proplot_cycles``,
        ``proplot_colors``, and ``proplot_fonts`` in this directory and in parent
        directories. "Hidden" folders with a leading dot are also recognized. Files
        in these directories are automatically loaded when proplot is imported.

        See also
        --------
        Configurator.user_folder
        Configurator.local_files
        """
        cdir = os.getcwd()
        paths = []
        if subfolder is None:
            subfolder = ("cmaps", "cycles", "colors", "fonts")
        if isinstance(subfolder, str):
            subfolder = (subfolder,)
        while cdir:  # i.e. not root
            for prefix in ("proplot", ".proplot"):
                for suffix in subfolder:
                    path = os.path.join(cdir, "_".join((prefix, suffix)))
                    if os.path.isdir(path):
                        paths.append(path)
            ndir = os.path.dirname(cdir)
            if ndir == cdir:  # root
                break
            cdir = ndir
        return paths[::-1]

    @staticmethod
    def _config_folder():
        """
        Get the XDG proplot folder.
        """
        home = os.path.expanduser("~")
        base = os.environ.get("XDG_CONFIG_HOME")
        if not base:
            base = os.path.join(home, ".config")
        if sys.platform.startswith(("linux", "freebsd")) and os.path.isdir(base):
            configdir = os.path.join(base, "proplot")
        else:
            configdir = os.path.join(home, ".proplot")
        return configdir

    @staticmethod
    def user_file():
        """
        Return location of the default proplotrc file. On Linux, this is either
        ``$XDG_CONFIG_HOME/proplot/proplotrc`` or ``~/.config/proplot/proplotrc``
        if the `XDG directory <https://wiki.archlinux.org/title/XDG_Base_Directory>`__
        is unset. On other operating systems, this is ``~/.proplot/proplotrc``. The
        location ``~/.proplotrc`` or ``~/.proplot/proplotrc`` is always returned if the
        file exists, regardless of the operating system. If multiple valid locations
        are found, a warning is raised.

        See also
        --------
        Configurator.user_folder
        Configurator.local_files
        """
        # Support both loose files and files inside .proplot
        file = os.path.join(Configurator.user_folder(), "proplotrc")
        universal = os.path.join(os.path.expanduser("~"), ".proplotrc")
        if os.path.isfile(universal):
            if file != universal and os.path.isfile(file):
                warnings._warn_proplot(
                    "Found conflicting default user proplotrc files at "
                    f"{universal!r} and {file!r}. Ignoring the second one."
                )
            file = universal
        return file

    @staticmethod
    def user_folder(subfolder=None):
        """
        Return location of the default proplot folder. On Linux, this
        is either ``$XDG_CONFIG_HOME/proplot`` or ``~/.config/proplot``
        if the `XDG directory <https://wiki.archlinux.org/title/XDG_Base_Directory>`__
        is unset. On other operating systems, this is ``~/.proplot``. The location
        ``~/.proplot`` is always returned if the folder exists, regardless of the
        operating system. If multiple valid locations are found, a warning is raised.

        See also
        --------
        Configurator.user_file
        Configurator.local_folders
        """
        # Try the XDG standard location
        # NOTE: This is borrowed from matplotlib.get_configdir
        home = os.path.expanduser("~")
        universal = folder = os.path.join(home, ".proplot")
        if sys.platform.startswith(("linux", "freebsd")):
            xdg = os.environ.get("XDG_CONFIG_HOME")
            xdg = xdg or os.path.join(home, ".config")
            folder = os.path.join(xdg, "proplot")
        # Fallback to the loose ~/.proplot if it is present
        # NOTE: This is critical or we might ignore previously stored settings!
        if os.path.isdir(universal):
            if folder != universal and os.path.isdir(folder):
                warnings._warn_proplot(
                    "Found conflicting default user proplot folders at "
                    f"{universal!r} and {folder!r}. Ignoring the second one."
                )
            folder = universal
        # Return the folder
        if subfolder:
            folder = os.path.join(folder, subfolder)
        return folder

    def context(self, *args, mode=0, file=None, **kwargs):
        """
        Temporarily modify the rc settings in a "with as" block.

        Parameters
        ----------
        *args
            Dictionaries of `rc` keys and values.
        file : path-like, optional
            Filename from which settings should be loaded.
        **kwargs
            `rc` names and values passed as keyword arguments.
            If the name has dots, simply omit them.

        Other parameters
        ----------------
        mode : {0, 1, 2}, optional
            The context mode. Dictates the behavior of `~Configurator.find`,
            `~Configurator.fill`, and `~Configurator.category` within a
            "with as" block when called with ``context=True``.

            The options are as follows:

            * ``mode=0``: Matplotlib's `rc_matplotlib` settings
              and proplot's `rc_proplot` settings are all returned,
              whether or not they are local to the "with as" block.
            * ``mode=1``: Matplotlib's `rc_matplotlib` settings are only
              returned if they are local to the "with as" block. For example,
              if :rcraw:`axes.titlesize` was passed to `~Configurator.context`,
              then ``pplt.rc.find('axes.titlesize', context=True)`` will return
              this value, but ``pplt.rc.find('axes.titleweight', context=True)`` will
              return ``None``. This is used internally when instantiating axes.
            * ``mode=2``: Matplotlib's `rc_matplotlib` settings and proplot's
              `rc_proplot` settings are only returned if they are local to the
              "with as" block. This is used internally when formatting axes.

        Note
        ----
        Context "modes" are primarily used internally but may also be useful for power
        users. Mode ``1`` is used when `~proplot.axes.Axes.format` is called during
        axes instantiation, and mode ``2`` is used when `~proplot.axes.Axes.format`
        is manually called by users. The latter prevents successive calls to
        `~proplot.axes.Axes.format` from constantly looking up and re-applying
        unchanged settings and significantly increasing the runtime.

        Example
        -------
        The below applies settings to axes in a specific figure using
        `~Configurator.context`.

        >>> import proplot as pplt
        >>> with pplt.rc.context(ticklen=5, metalinewidth=2):
        >>>     fig, ax = pplt.subplots()
        >>>     ax.plot(data)

        The below applies settings to a specific axes using
        `~proplot.axes.Axes.format`, which uses `~Configurator.context`
        internally.

        >>> import proplot as pplt
        >>> fig, ax = pplt.subplots()
        >>> ax.format(ticklen=5, metalinewidth=2)
        """
        # Add input dictionaries
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError(f"Non-dictionary argument {arg!r}.")
            kwargs.update(arg)

        # Add settings from file
        if file is not None:
            kw = self._load_file(file)
            kw = {key: value for key, value in kw.items() if key not in kwargs}
            kwargs.update(kw)

        # Activate context object
        if mode not in range(3):
            raise ValueError(f"Invalid mode {mode!r}.")
        cls = namedtuple("RcContext", ("mode", "kwargs", "rc_new", "rc_old"))
        context = cls(mode=mode, kwargs=kwargs, rc_new={}, rc_old={})
        self._context.append(context)
        return self

    def category(self, cat, *, trimcat=True, context=False):
        """
        Return a dictionary of settings beginning with the substring ``cat + '.'``.
        Optionally limit the search to the context level.

        Parameters
        ----------
        cat : str, optional
            The `rc` setting category.
        trimcat : bool, default: True
            Whether to trim ``cat`` from the key names in the output dictionary.
        context : bool, default: False
            If ``True``, then settings not found in the context dictionaries
            are omitted from the output dictionary. See `~Configurator.context`.

        See also
        --------
        Configurator.find
        Configurator.fill
        """
        kw = {}
        if cat not in rcsetup._rc_categories:
            raise ValueError(
                f"Invalid rc category {cat!r}. Valid categories are: "
                + ", ".join(map(repr, rcsetup._rc_categories))
                + "."
            )
        for key in self:
            if not re.match(rf"\A{cat}\.[^.]+\Z", key):
                continue
            value = self._get_item_context(key, None if context else 0)
            if value is None:
                continue
            if trimcat:
                key = re.sub(rf"\A{cat}\.", "", key)
            kw[key] = value
        return kw

    def fill(self, props, *, context=False):
        """
        Return a dictionary filled with settings whose names match the string values
        in the input dictionary. Optionally limit the search to the context level.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are setting names -- for example
            ``rc.fill({'edgecolor': 'axes.edgecolor', 'facecolor': 'axes.facecolor'})``.
        context : bool, default: False
            If ``True``, then settings not found in the context dictionaries
            are omitted from the output dictionary. See `~Configurator.context`.

        See also
        --------
        Configurator.category
        Configurator.find
        """
        kw = {}
        for key, value in props.items():
            item = self._get_item_context(value, None if context else 0)
            if item is not None:
                kw[key] = item
        return kw

    def find(self, key, *, context=False):
        """
        Return a single setting. Optionally limit the search to the context level.

        Parameters
        ----------
        key : str
            The single setting name.
        context : bool, default: False
            If ``True``, then ``None`` is returned if the setting is not found
            in the context dictionaries. See `~Configurator.context`.

        See also
        --------
        Configurator.category
        Configurator.fill
        """
        return self._get_item_context(key, None if context else 0)

    def update(self, *args, **kwargs):
        """
        Update several settings at once.

        Parameters
        ----------
        *args : str or dict-like, optional
            A dictionary containing `rc` keys and values. You can also pass
            a "category" name as the first argument, in which case all
            settings are prepended with ``'category.'``. For example,
            ``rc.update('axes', labelsize=20, titlesize=20)`` changes the
            :rcraw:`axes.labelsize` and :rcraw:`axes.titlesize` settings.
        **kwargs
            `rc` keys and values passed as keyword arguments.
            If the name has dots, simply omit them.

        See also
        --------
        Configurator.category
        Configurator.fill
        """
        prefix, kw = "", {}
        if not args:
            pass
        elif len(args) == 1 and isinstance(args[0], str):
            prefix = args[0]
        elif len(args) == 1 and isinstance(args[0], dict):
            kw = args[0]
        elif len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], dict):
            prefix, kw = args
        else:
            raise ValueError(
                f"Invalid arguments {args!r}. Usage is either "
                "rc.update(dict), rc.update(kwy=value, ...), "
                "rc.update(category, dict), or rc.update(category, key=value, ...)."
            )
        prefix = prefix and prefix + "."
        kw.update(kwargs)
        for key, value in kw.items():
            self.__setitem__(prefix + key, value)

    @docstring._snippet_manager
    def reset(self, local=True, user=True, default=True, **kwargs):
        """
        Reset the configurator to its initial state.

        Parameters
        ----------
        %(rc.params)s
        """
        self._init(local=local, user=user, default=default, **kwargs)

    def _load_file(self, path):
        """
        Return dictionaries of proplot and matplotlib settings loaded from the file.
        """
        # WARNING: Critical to not yet apply _get_item_dicts() syncing or else we
        # can overwrite input settings (e.g. label.size followed by font.size).
        path = os.path.expanduser(path)
        added = set()
        rcdict = {}
        with open(path, "r") as fh:
            for idx, line in enumerate(fh):
                # Strip comments
                message = f"line #{idx + 1} in file {path!r}"
                stripped = line.split("#", 1)[0].strip()
                if not stripped:
                    pass  # no warning
                    continue
                # Parse the pair
                pair = stripped.split(":", 1)
                if len(pair) != 2:
                    warnings._warn_proplot(f'Illegal {message}:\n{line}"')
                    continue
                # Detect duplicates
                key, value = map(str.strip, pair)
                if key in added:
                    warnings._warn_proplot(f"Duplicate rc key {key!r} on {message}.")
                added.add(key)
                # Get child dictionaries. Careful to have informative messages
                with warnings.catch_warnings():
                    warnings.simplefilter("error", warnings.ProplotWarning)
                    try:
                        key, value = self._validate_key(key, value)
                        value = self._validate_value(key, value)
                    except KeyError:
                        warnings.simplefilter("default", warnings.ProplotWarning)
                        warnings._warn_proplot(f"Invalid rc key {key!r} on {message}.")
                        continue
                    except ValueError as err:
                        warnings.simplefilter("default", warnings.ProplotWarning)
                        warnings._warn_proplot(
                            f"Invalid rc value {value!r} for key {key!r} on {message}: {err}"
                        )  # noqa: E501
                        continue
                    except warnings.ProplotWarning as err:
                        warnings.simplefilter("default", warnings.ProplotWarning)
                        warnings._warn_proplot(
                            f"Outdated rc key {key!r} on {message}: {err}"
                        )  # noqa: E501
                        warnings.simplefilter("ignore", warnings.ProplotWarning)
                        key, value = self._validate_key(key, value)
                        value = self._validate_value(key, value)
                # Update the settings
                rcdict[key] = value

        return rcdict

    def load(self, path):
        """
        Load settings from the specified file.

        Parameters
        ----------
        path : path-like
            The file path.

        See also
        --------
        Configurator.save
        """
        rcdict = self._load_file(path)
        for key, value in rcdict.items():
            self.__setitem__(key, value)

    @staticmethod
    def _save_rst(path):
        """
        Create an RST table file. Used for online docs.
        """
        string = rcsetup._rst_table()
        with open(path, "w") as fh:
            fh.write(string)

    @staticmethod
    def _save_yaml(path, user_dict=None, *, comment=False, description=False):
        """
        Create a YAML file. Used for online docs and default and user-generated
        proplotrc files. Extra settings can be passed with the input dictionary.
        """
        user_table = ()
        if user_dict:  # add always-uncommented user settings
            user_table = rcsetup._yaml_table(user_dict, comment=False)
            user_table = ("# Changed settings", user_table, "")
        proplot_dict = (
            rcsetup._rc_proplot_table if description else rcsetup._rc_proplot_default
        )  # noqa: E501
        proplot_table = rcsetup._yaml_table(
            proplot_dict, comment=comment, description=description
        )  # noqa: E501
        proplot_table = ("# Proplot settings", proplot_table, "")
        matplotlib_dict = rcsetup._rc_matplotlib_default
        matplotlib_table = rcsetup._yaml_table(matplotlib_dict, comment=comment)
        matplotlib_table = ("# Matplotlib settings", matplotlib_table)
        parts = (
            "#--------------------------------------------------------------------",
            "# Use this file to change the default proplot and matplotlib settings.",
            "# The syntax is identical to matplotlibrc syntax. For details see:",
            "# https://proplot.readthedocs.io/en/latest/configuration.html",
            "# https://matplotlib.org/stable/tutorials/introductory/customizing.html",
            "#--------------------------------------------------------------------",
            *user_table,  # empty if nothing passed
            *proplot_table,
            *matplotlib_table,
        )
        with open(path, "w") as fh:
            fh.write("\n".join(parts))

    def save(self, path=None, user=True, comment=None, backup=True, description=False):
        """
        Save the current settings to a ``proplotrc`` file. This writes
        the default values commented out plus the values that *differ*
        from the defaults at the top of the file.

        Parameters
        ----------
        path : path-like, default: 'proplotrc'
            The file name and/or directory. The default file name is ``proplotrc``
            and the default directory is the current directory.
        user : bool, default: True
            If ``True`` then settings that have been `~Configurator.changed` from
            the proplot defaults are shown uncommented at the top of the file.
        backup : bool, default: True
            Whether to "backup" an existing file by renaming with the suffix ``.bak``
            or overwrite an existing file.
        comment : bool, optional
            Whether to comment out the default settings. If not passed
            this takes the same value as `user`.
        description : bool, default: False
            Whether to include descriptions of each setting (as seen in the
            :ref:`user guide table <ug_rctable>`) as comments.

        See also
        --------
        Configurator.load
        Configurator.changed
        """
        path = os.path.expanduser(path or ".")
        if os.path.isdir(path):  # includes ''
            path = os.path.join(path, "proplotrc")
        if os.path.isfile(path) and backup:
            backup = path + ".bak"
            os.rename(path, backup)
            warnings._warn_proplot(f"Existing file {path!r} was moved to {backup!r}.")
        comment = _not_none(comment, user)
        user_dict = self.changed if user else None
        self._save_yaml(path, user_dict, comment=comment, description=description)

    @property
    def _context_mode(self):
        """
        Return the highest (least permissive) context mode.
        """
        return max((context.mode for context in self._context), default=0)

    @property
    def changed(self):
        """
        A dictionary of settings that have changed from the proplot defaults.

        See also
        --------
        Configurator.save
        """
        # Carefully detect changed settings
        rcdict = {}
        for key, value in self.items():
            default = rcsetup._get_default_param(key)
            if (
                isinstance(value, Real)
                and isinstance(default, Real)
                and np.isclose(value, default)
            ):  # noqa: E501
                pass
            elif value == default:
                pass
            else:
                rcdict[key] = value
        # Ignore non-style-related settings. See mstyle.STYLE_BLACKLIST
        # TODO: For now not sure how to detect if prop cycle changed since
        # we cannot load it from _cmap_database in rcsetup.
        rcdict.pop("interactive", None)  # changed by backend
        rcdict.pop("axes.prop_cycle", None)
        return _filter_style_dict(rcdict, warn=False)

    # Renamed methods
    load_file = warnings._rename_objs("0.8.0", load_file=load)


# Initialize locations
_init_user_folders()
_init_user_file()

#: A dictionary-like container of matplotlib settings. Assignments are
#: validated and restricted to recognized setting names.
rc_matplotlib = mpl.rcParams  # PEP8 4 lyfe

#: A dictionary-like container of proplot settings. Assignments are
#: validated and restricted to recognized setting names.
rc_proplot = rcsetup._rc_proplot_default.copy()  # a validated rcParams-style dict

#: Instance of `Configurator`. This controls both `rc_matplotlib` and `rc_proplot`
#: settings. See the :ref:`configuration guide <ug_config>` for details.
rc = Configurator(skip_cycle=True)

# Deprecated
RcConfigurator = warnings._rename_objs(
    "0.8.0",
    RcConfigurator=Configurator,
)
inline_backend_fmt = warnings._rename_objs(
    "0.6.0", inline_backend_fmt=config_inline_backend
)
