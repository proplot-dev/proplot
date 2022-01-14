#!/usr/bin/env python3
"""
Custom :rc: and :rcraw: roles for rc settings.
"""
from docutils import nodes
from os.path import sep
from proplot.internals import rcsetup
from matplotlib import rcParams


def _get_nodes(rawtext, text, inliner):
    """
    Return the literal node.
    """
    source = inliner.document.attributes['source'].replace(sep, '/')
    relsource = source.split('/docs/', 1)
    if len(relsource) == 1:
        return []
    if text in rcParams:
        refuri = 'https://matplotlib.org/stable/tutorials/introductory/customizing.html'
        refuri = f'{refuri}?highlight={text}#the-matplotlibrc-file'
    else:
        path = '../' * relsource[1].count('/') + 'en/stable'
        refuri = f'{path}/configuration.html?highlight={text}#table-of-settings'
    node = nodes.Text(f'rc[{text!r}]' if '.' in text else f'rc.{text}')
    ref = nodes.reference(rawtext, node, refuri=refuri)
    return [nodes.literal('', '', ref)]


def rc_name_role(name, rawtext, text, lineno, inliner, options={}, content=[]):  # noqa: U100, E501
    """
    The :rcname: role. Includes a link to the setting.
    """
    list_ = _get_nodes(rawtext, text, inliner)
    return list_, []


def rc_role(name, rawtext, text, lineno, inliner, options={}, content=[]):  # noqa: U100
    """
    The :rc: role. Includes a link to the setting and its default value.
    """
    list_ = _get_nodes(rawtext, text, inliner)
    try:
        default = rcsetup._get_default_param(text)
    except KeyError:
        pass
    else:
        list_.append(nodes.Text(' = '))
        list_.append(nodes.literal('', '', nodes.Text(repr(default))))
    return nodes, []


def setup(app):
    """
    Set up the roles.
    """
    app.add_role('rc', rc_role)
    app.add_role('rcname', rc_name_role)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
