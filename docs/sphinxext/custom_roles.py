#!/usr/bin/env python3
"""
Custom :rc: and :rcraw: roles for rc settings.
"""
import os

from docutils import nodes
from matplotlib import rcParams

from proplot.internals import rcsetup


def _node_list(rawtext, text, inliner):
    """
    Return a singleton node list or an empty list if source is unknown.
    """
    source = inliner.document.attributes["source"].replace(os.path.sep, "/")
    relsource = source.split("/docs/", 1)
    if len(relsource) == 1:
        return []
    if text in rcParams:
        refuri = "https://matplotlib.org/stable/tutorials/introductory/customizing.html"
        refuri = f"{refuri}?highlight={text}#the-matplotlibrc-file"
    else:
        path = "../" * relsource[1].count("/") + "en/stable"
        refuri = f"{path}/configuration.html?highlight={text}#table-of-settings"
    node = nodes.Text(f"rc[{text!r}]" if "." in text else f"rc.{text}")
    ref = nodes.reference(rawtext, node, refuri=refuri)
    return [nodes.literal("", "", ref)]


def rc_raw_role(
    name, rawtext, text, lineno, inliner, options={}, content=[]
):  # noqa: U100, E501
    """
    The :rcraw: role. Includes a link to the setting.
    """
    node_list = _node_list(rawtext, text, inliner)
    return node_list, []


def rc_role(name, rawtext, text, lineno, inliner, options={}, content=[]):  # noqa: U100
    """
    The :rc: role. Includes a link to the setting and its default value.
    """
    node_list = _node_list(rawtext, text, inliner)
    try:
        default = rcsetup._get_default_param(text)
    except KeyError:
        pass
    else:
        node_list.append(nodes.Text(" = "))
        node_list.append(nodes.literal("", "", nodes.Text(repr(default))))
    return node_list, []


def setup(app):
    """
    Set up the roles.
    """
    app.add_role("rc", rc_role)
    app.add_role("rcraw", rc_raw_role)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
