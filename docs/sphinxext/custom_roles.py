from docutils import nodes
from os.path import sep
from proplot.internals import rcsetup


def get_nodes(rawtext, text, inliner):
    repr_ = f"rc['{text}']" if '.' in text else f'rc.{text}'
    section = 'rc-added' if '.' in text else 'rc-quick'  # reST labels for sections
    rendered = nodes.Text(repr_)
    source = inliner.document.attributes['source'].replace(sep, '/')
    relsource = source.split('/docs/', 1)
    if len(relsource) == 1:
        return []
    levels = relsource[1].count('/')  # distance to 'docs' folder
    refuri = '../' * levels + f'en/latest/configuration.html?highlight={text}#{section}'
    ref = nodes.reference(rawtext, rendered, refuri=refuri)
    return [nodes.literal('', '', ref)]


def rc_role(name, rawtext, text, lineno, inliner, options={}, content=[]):  # noqa: U100
    node_list = get_nodes(rawtext, text, inliner)
    try:
        default = rcsetup._get_default_param(text)
    except KeyError:
        pass
    else:
        node_list.append(nodes.Text(' = '))
        node_list.append(nodes.literal('', '', nodes.Text(repr(default))))
    return node_list, []


def rc_role_raw(name, rawtext, text, lineno, inliner, options={}, content=[]):  # noqa: U100, E501
    return get_nodes(rawtext, text, inliner), []


def setup(app):
    app.add_role('rc', rc_role)
    app.add_role('rcraw', rc_role_raw)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
