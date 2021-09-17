from docutils import nodes
from os.path import sep
from proplot.internals import rcsetup
from matplotlib import rcParams


def get_nodes(rawtext, text, inliner):
    repr_ = f"rc['{text}']" if '.' in text else f'rc.{text}'
    rendered = nodes.Text(repr_)
    source = inliner.document.attributes['source'].replace(sep, '/')
    relsource = source.split('/docs/', 1)
    if len(relsource) == 1:
        return []
    if text in rcParams:
        refuri = 'https://matplotlib.org/stable/tutorials/introductory/customizing.html'
        refuri = f'{refuri}?highlight={text}#the-matplotlibrc-file'
    else:
        levels = relsource[1].count('/')  # distance to base URL
        refuri = '../' * levels + 'en/stable/configuration.html'
        refuri += f'?highlight={text}#rc-proplot'
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
