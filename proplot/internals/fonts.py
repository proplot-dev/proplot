#!/usr/bin/env python3
"""
Overrides related to math fonts
"""
import matplotlib as mpl
from matplotlib.font_manager import findfont, ttfFontProperty
from matplotlib.mathtext import MathTextParser

from . import warnings

try:  # newer versions
    from matplotlib._mathtext import UnicodeFonts
except ImportError:  # older versions
    from matplotlib.mathtext import UnicodeFonts


class _UnicodeFonts(UnicodeFonts):
    """
    A simple `~matplotlib._mathtext.UnicodeFonts` subclass that
    interprets ``rc['mathtext.default'] != 'regular'`` in the presence of
    ``rc['mathtext.fontset'] == 'custom'`` as possibly modifying the active font.

    Works by permitting the ``rc['mathtext.rm']``, ``rc['mathtext.it']``,
    etc. settings to have the dummy value ``'current'`` instead of a valid family
    name, e.g. ``rc['mathtext.it'] == 'current:italic'`` (permitted through an
    override of the `~matplotlib.rcsetup.validate_font_properties` validator).
    When this dummy value is detected then the font properties passed to
    `~matplotlib._mathtext.TrueTypeFont` are taken by replacing ``'current'``
    in the "math" fontset with the active font name.
    """
    def __init__(self, *args, **kwargs):
        # Initialize font
        # NOTE: Could also capture the 'default_font_prop' passed as positional
        # argument but want to guard against keyword changes. This entire API is
        # private and it is easier to do graceful fallback with _fonts dictionary.
        ctx = {}  # rc context
        regular = {}  # styles
        for texfont in ('cal', 'rm', 'tt', 'it', 'bf', 'sf'):
            key = 'mathtext.' + texfont
            prop = mpl.rcParams[key]
            if prop.startswith('regular'):
                ctx[key] = prop.replace('regular', 'sans', 1)
                regular[texfont] = prop
        with mpl.rc_context(ctx):
            super().__init__(*args, **kwargs)
        # Apply current font replacements
        if (
            not hasattr(self, 'fontmap')
            or not hasattr(self, '_fonts')
            or 'regular' not in self._fonts
        ):
            warnings._warn_proplot('Failed to update math text parser.')
        else:
            font = self._fonts['regular']  # an ft2font.FT2Font instance
            font = ttfFontProperty(font)
            for texfont, prop in regular.items():
                prop = prop.replace('regular', font.name)
                self.fontmap[texfont] = findfont(prop, fallback_to_default=False)


# Replace the parser
try:
    mapping = MathTextParser._font_type_mapping
    if mapping['custom'] is UnicodeFonts:
        mapping['custom'] = _UnicodeFonts
except (KeyError, AttributeError):
    warnings._warn_proplot('Failed to update math text parser.')
