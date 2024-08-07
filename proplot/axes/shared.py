#!/usr/bin/env python
"""
An axes used to jointly format Cartesian and polar axes.
"""
# NOTE: We could define these in base.py but idea is projection-specific formatters
# should never be defined on the base class. Might add to this class later anyway.
import numpy as np

from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import _pop_kwargs
from ..utils import _fontsize_to_pt, _not_none, units


class _SharedAxes(object):
    """
    Mix-in class with methods shared between `~proplot.axes.CartesianAxes`
    and `~proplot.axes.PolarAxes`.
    """

    @staticmethod
    def _min_max_lim(key, min_=None, max_=None, lim=None):
        """
        Translate and standardize minimum, maximum, and limit keyword arguments.
        """
        if lim is None:
            lim = (None, None)
        if not np.iterable(lim) or not len(lim) == 2:
            raise ValueError(f"Invalid {key}{lim!r}. Must be 2-tuple of values.")
        min_ = _not_none(**{f"{key}min": min_, f"{key}lim_0": lim[0]})
        max_ = _not_none(**{f"{key}max": max_, f"{key}lim_1": lim[1]})
        return min_, max_

    def _update_background(self, x=None, tickwidth=None, tickwidthratio=None, **kwargs):
        """
        Update the background patch and spines.
        """
        # Update the background patch
        kw_face, kw_edge = rc._get_background_props(**kwargs)
        self.patch.update(kw_face)
        if x is None:
            opts = self.spines
        elif x == "x":
            opts = ("bottom", "top", "inner", "polar")
        else:
            opts = ("left", "right", "start", "end")
        for opt in opts:
            self.spines.get(opt, {}).update(kw_edge)

        # Update the tick colors
        axis = "both" if x is None else x
        x = _not_none(x, "x")
        obj = getattr(self, x + "axis")
        edgecolor = kw_edge.get("edgecolor", None)
        if edgecolor is not None:
            self.tick_params(axis=axis, which="both", color=edgecolor)

        # Update the tick widths
        # NOTE: Only use 'linewidth' if it was explicitly passed. Do not
        # include 'linewidth' inferred from rc['axes.linewidth'] setting.
        kwmajor = getattr(obj, "_major_tick_kw", {})  # graceful fallback if API changes
        kwminor = getattr(obj, "_minor_tick_kw", {})
        if "linewidth" in kwargs:
            tickwidth = _not_none(tickwidth, kwargs["linewidth"])
        tickwidth = _not_none(tickwidth, rc.find("tick.width", context=True))
        tickwidthratio = _not_none(
            tickwidthratio, rc.find("tick.widthratio", context=True)
        )  # noqa: E501
        tickwidth_prev = kwmajor.get("width", rc[x + "tick.major.width"])
        if tickwidth_prev == 0:
            tickwidthratio_prev = rc["tick.widthratio"]  # no other way of knowing
        else:
            tickwidthratio_prev = (
                kwminor.get("width", rc[x + "tick.minor.width"]) / tickwidth_prev
            )  # noqa: E501
        for which in ("major", "minor"):
            kwticks = {}
            if tickwidth is not None or tickwidthratio is not None:
                tickwidth = _not_none(tickwidth, tickwidth_prev)
                kwticks["width"] = tickwidth = units(tickwidth, "pt")
                if tickwidth == 0:  # avoid unnecessary padding
                    kwticks["size"] = 0
                elif which == "minor":
                    tickwidthratio = _not_none(tickwidthratio, tickwidthratio_prev)
                    kwticks["width"] *= tickwidthratio
            self.tick_params(axis=axis, which=which, **kwticks)

    def _update_ticks(
        self,
        x,
        *,
        grid=None,
        gridminor=None,
        gridpad=None,
        gridcolor=None,
        ticklen=None,
        ticklenratio=None,
        tickdir=None,
        tickcolor=None,
        labeldir=None,
        labelpad=None,
        labelcolor=None,
        labelsize=None,
        labelweight=None,
    ):
        """
        Update the gridlines and labels. Set `gridpad` to ``True`` to use grid padding.
        """
        # Filter out text properties
        axis = "both" if x is None else x
        kwtext = rc._get_ticklabel_props(axis)
        kwtext_extra = _pop_kwargs(kwtext, "weight", "family")
        kwtext = {"label" + key: value for key, value in kwtext.items()}
        if labelcolor is not None:
            kwtext["labelcolor"] = labelcolor
        if labelsize is not None:
            kwtext["labelsize"] = labelsize
        if labelweight is not None:
            kwtext_extra["weight"] = labelweight

        # Apply tick settings with tick_params when possible
        x = _not_none(x, "x")
        obj = getattr(self, x + "axis")
        kwmajor = getattr(obj, "_major_tick_kw", {})  # graceful fallback if API changes
        kwminor = getattr(obj, "_minor_tick_kw", {})
        ticklen_prev = kwmajor.get("size", rc[x + "tick.major.size"])
        if ticklen_prev == 0:
            ticklenratio_prev = rc["tick.lenratio"]  # no other way of knowing
        else:
            ticklenratio_prev = (
                kwminor.get("size", rc[x + "tick.minor.size"]) / ticklen_prev
            )  # noqa: E501
        for b, which in zip((grid, gridminor), ("major", "minor")):
            # Tick properties
            # NOTE: Must make 'tickcolor' overwrite 'labelcolor' or else 'color'
            # passed to __init__ will not apply correctly. Annoying but unavoidable
            kwticks = rc._get_tickline_props(axis, which=which)
            if labelpad is not None:
                kwticks["pad"] = labelpad
            if tickcolor is not None:
                kwticks["color"] = tickcolor
            if ticklen is not None or ticklenratio is not None:
                ticklen = _not_none(ticklen, ticklen_prev)
                kwticks["size"] = ticklen = units(ticklen, "pt")
                if ticklen > 0 and which == "minor":
                    ticklenratio = _not_none(ticklenratio, ticklenratio_prev)
                    kwticks["size"] *= ticklenratio
            if gridpad:  # use grid.labelpad instead of tick.labelpad
                kwticks.pop("pad", None)
                pad = rc.find("grid.labelpad", context=True)
                if pad is not None:
                    kwticks["pad"] = units(pad, "pt")

            # Tick direction properties
            # NOTE: These have no x and y-specific versions but apply here anyway
            if labeldir == "in":  # put tick labels inside the plot
                tickdir = "in"
                kwticks.setdefault(
                    "pad",
                    -rc[f"{axis}tick.major.size"]
                    - _not_none(labelpad, rc[f"{axis}tick.major.pad"])
                    - _fontsize_to_pt(rc[f"{axis}tick.labelsize"]),
                )
            if tickdir is not None:
                kwticks["direction"] = tickdir

            # Gridline properties
            # NOTE: Internally ax.grid() passes gridOn to ax.tick_params() but this
            # is undocumented and might have weird side effects. Just use ax.grid()
            b = rc._get_gridline_bool(b, axis=axis, which=which)
            if b is not None:
                self.grid(b, axis=axis, which=which)
            kwlines = rc._get_gridline_props(which=which)
            if "axisbelow" in kwlines:
                self.set_axisbelow(kwlines.pop("axisbelow"))
            if gridcolor is not None:
                kwlines["grid_color"] = gridcolor

            # Apply tick and gridline properties
            kwticks.pop("ndivs", None)  # not in mpl
            self.tick_params(axis=axis, which=which, **kwticks, **kwlines, **kwtext)

        # Apply settings that can't be controlled with tick_params
        if kwtext_extra:
            for lab in obj.get_ticklabels():
                lab.update(kwtext_extra)
